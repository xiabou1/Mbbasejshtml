// mandelbulbWorker.js

// Funciones matemáticas básicas
function normalize(v) {
    let l = Math.hypot(v[0], v[1], v[2]);
    if (l === 0) return [0,0,0];
    return [v[0]/l, v[1]/l, v[2]/l];
}

function dot(a, b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

function cross(a,b) {
    return [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ];
}

function add(a,b) {
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2]];
}

function sub(a,b) {
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
}

function mul(v, s) {
    return [v[0]*s, v[1]*s, v[2]*s];
}

// Funciones de Mandelbulb y Raymarching
let maxMandelbulbIterations_worker;
let maxRaymarchingSteps_worker;

function mandelbulb(p) {
    let z = [...p];
    let dr = 1.0;
    let r = 0.0;
    const Power = 8;
    for (let i = 0; i < maxMandelbulbIterations_worker; i++) {
        r = Math.hypot(z[0], z[1], z[2]);
        if (r > 2.0) break;

        let theta = Math.acos(z[2]/r);
        let phi = Math.atan2(z[1], z[0]);

        dr = Math.pow(r, Power - 1.0) * Power * dr + 1.0;

        let zr = Math.pow(r, Power);
        theta *= Power;
        phi *= Power;

        z[0] = zr * Math.sin(theta) * Math.cos(phi) + p[0];
        z[1] = zr * Math.sin(theta) * Math.sin(phi) + p[1];
        z[2] = zr * Math.cos(theta) + p[2];
    }
    return 0.5 * Math.log(r) * r / dr;
}

function traceRay(startPos, direction, maxSteps, minDistThreshold, fractalScale) {
    let p = [...startPos];
    let hit = false;
    let distance = 0;

    for (let i = 0; i < maxSteps; i++) {
        let d = mandelbulb(p);

        d *= fractalScale;

        if (d < minDistThreshold) {
            hit = true;
            break;
        }
        if (d > 4.0 * fractalScale) break; // Evita que los rayos viajen demasiado lejos
        p = add(p, mul(direction, d));
        distance += d;
    }
    return { hit: hit, point: p, distance: distance };
}

function calculateNormal(p, currentMinDist) {
    const epsilon = 0.0001 * currentMinDist;
    const dx = mandelbulb([p[0] + epsilon, p[1], p[2]]) - mandelbulb([p[0] - epsilon, p[1], p[2]]);
    const dy = mandelbulb([p[0], p[1] + epsilon, p[2]]) - mandelbulb([p[0], p[1] - epsilon, p[2]]);
    const dz = mandelbulb([p[0], p[1], p[2] + epsilon]) - mandelbulb([p[0], p[1] - epsilon, p[2]]);
    return normalize([dx, dy, dz]);
}

// Event listener para recibir mensajes del hilo principal
self.onmessage = function(e) {
    const {
        command,
        params,
        width,
        height,
        yStart,
        yEnd
    } = e.data;

    if (command === 'renderChunk') {
        const {
            cameraDistanceExp,
            fractalZoomExp,
            polar,
            azimuthal,
            panX,
            panY,
            minDistExp,
            pixelSize,
            renderMode,
            eyeSeparation,
            convergenceFactor,
            maxMandelbulbIterations,
            maxRaymarchingSteps,
            lastCentralHitPoint,
            lightingBalance
        } = params;

        maxMandelbulbIterations_worker = maxMandelbulbIterations;
        maxRaymarchingSteps_worker = maxRaymarchingSteps;

        const currentCameraDistance = 4.0 * Math.pow(10, -cameraDistanceExp);
        const currentFractalScale = Math.pow(1.5, -fractalZoomExp);
        const currentMinDist = Math.pow(1.5, -minDistExp);

        const chunkHeight = yEnd - yStart;
        const chunkImageData = new ImageData(width, chunkHeight);
        const data = chunkImageData.data;

        const vec3 = () => [0,0,0];
        const eyeBase = vec3();

        eyeBase[0] = currentCameraDistance * Math.sin(polar * Math.PI/180) * Math.cos(azimuthal * Math.PI/180);
        eyeBase[1] = currentCameraDistance * Math.cos(polar * Math.PI/180);
        eyeBase[2] = currentCameraDistance * Math.sin(polar * Math.PI/180) * Math.sin(azimuthal * Math.PI/180);
        
        let cameraPosition = add(eyeBase, [panX, panY, 0]);

        let cameraLookAtTarget = [0,0,0];

        let initialUp = [0, 1, 0];

        let forwardBase = normalize(sub(cameraLookAtTarget, cameraPosition));
        let rightBase = normalize(cross(forwardBase, initialUp));
        let trueUpBase = cross(rightBase, forwardBase);

        let fov = 1.0;

        const light1Dir = normalize([-0.5, 0.8, -0.7]);
        const light1Color = [1.0, 0.8, 0.6];
        const light2Dir = normalize([0.7, 0.5, 0.5]);
        const light2Color = [0.6, 0.8, 1.0];
        const ambientStrength = 0.1;
        const specularStrength = 0.5;
        const shininess = 32.0;

        // Calcular los pesos dinámicos basados en lightingBalance
        const stepsContributionNormalized = (lightingBalance - 1) / 8.0; 
        const stepsWeight = 0.9 - (0.8 * stepsContributionNormalized); 
        const phongWeight = 1.0 - stepsWeight;

        let renderedPixelsInChunk = 0;

        for (let y_relative = 0; y_relative < chunkHeight; y_relative += pixelSize) {
            const y_absolute = yStart + y_relative;
            for (let x = 0; x < width; x += pixelSize) {
                let currentEye;
                let currentCenterTarget;
                let currentForward;
                let currentRight;
                let currentTrueUp;

                if (renderMode === "stereoscopic") {
                    const halfWidth = width / 2;
                    let offsetDirection = (x < halfWidth) ? -1 : 1;

                    let eyeOffset = mul(rightBase, offsetDirection * eyeSeparation);
                    currentEye = add(cameraPosition, eyeOffset);

                    let baseTargetForConvergence = lastCentralHitPoint || [0,0,0]; 
                    let centerOffset = mul(rightBase, offsetDirection * eyeSeparation * convergenceFactor);
                    currentCenterTarget = add(baseTargetForConvergence, centerOffset);

                    currentForward = normalize(sub(currentCenterTarget, currentEye));
                    currentRight = normalize(cross(currentForward, trueUpBase));
                    let potentialTrueUp = cross(currentRight, currentForward);
                    currentTrueUp = dot(potentialTrueUp, initialUp) > 0 ? potentialTrueUp : mul(potentialTrueUp, -1);

                } else {
                    currentEye = [...cameraPosition];
                    currentCenterTarget = [0,0,0];
                    currentForward = normalize(sub(currentCenterTarget, currentEye));
                    currentRight = normalize(cross(currentForward, initialUp));
                    currentTrueUp = cross(currentRight, currentForward);
                }

                let u = (x - width/2) / height * fov;
                let v = (y_absolute - height/2) / height * fov;

                let dir = normalize(add(add(currentForward, mul(currentRight, u * currentFractalScale)), mul(currentTrueUp, v * currentFractalScale)));

                const pixelHitResult = traceRay(currentEye, dir, maxRaymarchingSteps, currentMinDist, currentFractalScale);
                const hit = pixelHitResult.hit;
                const p = pixelHitResult.point;
                
                // --- AJUSTES CLAVE AQUÍ PARA INTENSIFICAR EL EFECTO DE LOS PASOS ---
                const stepsNormalized = pixelHitResult.distance / (4.0 * currentFractalScale);
                
                const rawStepsIntensity = Math.pow(1.0 - Math.min(1.0, stepsNormalized), 4.0); 
                
                const stepsIntensityMultiplier = 2.0; 
                const stepsIntensity = rawStepsIntensity * stepsIntensityMultiplier;

                let r_final = 0, g_final = 0, b_final = 0;

                if (hit) {
                    const normal = calculateNormal(p, currentMinDist);
                    const viewDir = normalize(sub(currentEye, p));

                    let ambient1 = mul(light1Color, ambientStrength);
                    let diff1 = Math.max(dot(normal, light1Dir), 0.0);
                    let diffuse1 = mul(light1Color, diff1);
                    let reflectDir1 = sub(mul(normal, 2 * dot(normal, light1Dir)), light1Dir);
                    let spec1 = Math.pow(Math.max(dot(viewDir, reflectDir1), 0.0), shininess);
                    let specular1 = mul(light1Color, spec1 * specularStrength);

                    let ambient2 = mul(light2Color, ambientStrength);
                    let diff2 = Math.max(dot(normal, light2Dir), 0.0);
                    let diffuse2 = mul(light2Color, diff2);
                    let reflectDir2 = sub(mul(normal, 2 * dot(normal, light2Dir)), light2Dir);
                    let spec2 = Math.pow(Math.max(dot(viewDir, reflectDir2), 0.0), shininess);
                    let specular2 = mul(light2Color, spec2 * specularStrength);

                    let phongR = (ambient1[0] + diffuse1[0] + specular1[0]) + (ambient2[0] + diffuse2[0] + specular2[0]);
                    let phongG = (ambient1[1] + diffuse1[1] + specular1[1]) + (ambient2[1] + diffuse2[1] + specular2[1]);
                    let phongB = (ambient1[2] + diffuse1[2] + specular1[2]) + (ambient2[2] + diffuse2[2] + specular2[2]);

                    const stepsBaseColor = [0.2, 0.5, 0.8]; 
                    let stepsColorR = stepsBaseColor[0] * stepsIntensity;
                    let stepsColorG = stepsBaseColor[1] * stepsIntensity;
                    let stepsColorB = stepsBaseColor[2] * stepsIntensity;

                    // Aplicar el balance de iluminación
                    r_final = Math.floor(((phongR * phongWeight) + (stepsColorR * stepsWeight)) * 255);
                    g_final = Math.floor(((phongG * phongWeight) + (stepsColorG * stepsWeight)) * 255);
                    b_final = Math.floor(((phongB * phongWeight) + (stepsColorB * stepsWeight)) * 255);

                    r_final = Math.min(255, Math.max(0, r_final));
                    g_final = Math.min(255, Math.max(0, g_final));
                    b_final = Math.min(255, Math.max(0, b_final));

                } else {
                    r_final = 0;
                    g_final = 0;
                    b_final = 0;
                }

                // Escribir en el ImageData del CHUNK
                for (let dy = 0; dy < pixelSize; dy++) {
                    for (let dx = 0; dx < pixelSize; dx++) {
                        if (x + dx < width && y_relative + dy < chunkHeight) {
                            let idx = 4 * ((y_relative + dy) * width + (x + dx));
                            data[idx + 0] = r_final;
                            data[idx + 1] = g_final;
                            data[idx + 2] = b_final;
                            data[idx + 3] = 255;
                        }
                    }
                }
                renderedPixelsInChunk += (pixelSize * pixelSize);
            }
        }
        self.postMessage({
            command: 'chunkRendered',
            imageData: chunkImageData,
            yStart: yStart,
            yEnd: yEnd,
            renderedPixelsInChunk: renderedPixelsInChunk
        }, [chunkImageData.data.buffer]);
    }
};

