import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';
import { OutputPass } from 'three/addons/postprocessing/OutputPass.js';

// --- グローバル変数 ---
let scene, camera, renderer, controls, composer;
let starParticles; // 遠景用のパーティクル
let nearStarsGroup; // 近景用のスプライトをまとめるグループ
let nearSphereGroup; // 最近景用の球体をまとめるグループ
let flareTexture; // フレア画像のテクスチャ
let starDataCache; // CSVから読み込んだ星の生データ
let originalSizes; // パーティクルの元のサイズを保持

// --- ブルーム効果のパラメータ ---
const bloomParams = {
    strength: 0.0,  // ブルームの強さ
    radius: 0.2,    // ブルームが広がる半径（小さくしてシャープに）
    threshold: 1.0, // この輝度以上のピクセルにブルームを適用
};

// --- 星の描画パラメータ ---
const starRenderParams = {
    minBrightness: 0.5,  // シェーダー内の最小輝度係数
    maxBrightness: 15.0, // シェーダー内の最大輝度係数
    sizeMultiplier: 0.15, // 星の基本サイズ係数
    minSize: 0.01,        // 星の最小描画サイズ
    maxSize: 25.0,       // 星の最大描画サイズ
    sizeLogBase: 10.0,   // 等級からサイズを計算する際の対数の底
    sizeLogFactor: -0.1, // 等級からサイズを計算する際の係数
    brightnessLogBase: 2.0, // 等級から明るさを計算する際の対数の底
    brightnessReferenceMag: 6.0, // 明るさ計算の基準等級
    flareDistanceFactor: 1.0, // 近景フレアの距離に応じた基本サイズ係数
    flareMinScale: 0.5,        // 近景フレアの最小スケール
    flareMaxScale: 20.0,       // 近景フレアの最大スケール
    sphereBaseScale: 0.005,     // 最近景の球体の基本スケール
    sphereMinScale: 0.01,      // 最近景の球体の最小スケール
    sphereMaxScale: 1.0,       // 最近景の球体の最大スケール
    distanceFalloffScale: 1600.0, // 遠景パーティクルの距離減衰係数 (基準距離40の2乗)
};

// --- LOD (Level of Detail) パラメータ ---
const lodDistance = 20.0; // この距離より近づいたらスプライトに切り替える
const lodSphereDistance = 0.5; // この距離より内側は完全に球体
const lodTransitionRange = 1.0; // 球体とフレアのクロスフェード範囲 (lodSphereDistanceからの距離)

// --- ヘルパー関数: スペクトル型から色を取得 ---
const getStarColor = (spectralType) => {
    // CSVから読み込んだ際に前後に空白が入る可能性を考慮してtrim()する
    const type = spectralType.trim();
    const colorMap = {
        "O": "#9bb0ff", "B": "#aabfff", "A": "#cad8ff",
        "F": "#f8f7ff", "G": "#fff4ea", "K": "#ffd2a1",
        "M": "#ffb56c", "WD": "#f0f0f0", "ES": "#ffeea8",
    };
    return colorMap[type] || "#ffffff";
};

// --- 初期化関数 ---
function init() {
    // シーン
    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x000005); // 少し青みがかった黒

    // カメラ
    camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 2000);
    camera.position.set(0, 50, 150); // 初期カメラ位置

    // レンダラー
    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(window.innerWidth, window.innerHeight);
    document.body.appendChild(renderer.domElement);

    // カメラコントロール
    controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.05;
    controls.screenSpacePanning = false;
    controls.minDistance = 1;
    controls.maxDistance = 500;

    // 近景用の星を管理するグループ
    nearStarsGroup = new THREE.Group();
    scene.add(nearStarsGroup);

    // 最近景用の星を管理するグループ
    nearSphereGroup = new THREE.Group();
    scene.add(nearSphereGroup);

    // ポストプロセッシング（ブルーム効果）
    const renderScene = new RenderPass(scene, camera);

    // ブルーム効果: パラメータを一元管理
    const bloomPass = new UnrealBloomPass(
        new THREE.Vector2(window.innerWidth, window.innerHeight),
        bloomParams.strength, bloomParams.radius, bloomParams.threshold
    );
    const outputPass = new OutputPass();

    composer = new EffectComposer(renderer);
    composer.addPass(renderScene);
    composer.addPass(bloomPass);
    composer.addPass(outputPass);
    // ウィンドウリサイズへの対応
    window.addEventListener('resize', onWindowResize, false);

    // データを読み込んで星を生成
    loadStarData('3d_viewer_data.csv');
}

// --- ウィンドウリサイズ処理 ---
function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
    composer.setSize(window.innerWidth, window.innerHeight);
}

// --- アニメーションループ ---
function animate() {
    requestAnimationFrame(animate);
    controls.update();
    updateStarLOD(); // 星のLODを更新
    composer.render(); // EffectComposerでレンダリング
}

// --- CSVデータから星を読み込む ---
async function loadStarData(csvPath) {
    try {
        const response = await fetch(csvPath);
        if (!response.ok) {
            throw new Error(`HTTP error! status: ${response.status}`);
        }
        const csvText = await response.text();
        const stars = parseCSV(csvText);
        starDataCache = stars; // データをキャッシュしておく
        await createStars(starDataCache);
    } catch (error) {
        console.error("Error loading or parsing star data:", error);
        const info = document.getElementById('info');
        info.textContent = "Error: Could not load 3d_viewer_data.csv. Make sure it has been generated and exists in the '3d_viewer' directory.";
        info.style.color = "red";
    }
}

// --- LODヘルパー: プールから球体を取得・更新 ---
function updateLODSphere(poolIndex, star, starPosition, opacity, magnitudeFactor) {
    let sphere = nearSphereGroup.children[poolIndex];
    if (!sphere) {
        const geometry = new THREE.SphereGeometry(1, 16, 16);
        const material = new THREE.MeshBasicMaterial({ transparent: true });
        sphere = new THREE.Mesh(geometry, material);
        nearSphereGroup.add(sphere);
    }
    sphere.visible = true;
    sphere.material.opacity = opacity;
    sphere.material.color.set(getStarColor(star.spectralType));
    sphere.position.copy(starPosition);

    // スケールを計算 (距離に依存せず、星の等級にのみ依存)
    const scale = starRenderParams.sphereBaseScale * magnitudeFactor;
    const finalScale = THREE.MathUtils.clamp(scale, starRenderParams.sphereMinScale, starRenderParams.sphereMaxScale);
    sphere.scale.set(finalScale, finalScale, finalScale);
}

// --- LODヘルパー: プールからスプライトを取得・更新 ---
function updateLODSprite(poolIndex, star, starPosition, opacity, magnitudeFactor, distance) {
    let sprite = nearStarsGroup.children[poolIndex];
    if (!sprite) {
        sprite = new THREE.Sprite(new THREE.SpriteMaterial({
            map: flareTexture,
            blending: THREE.AdditiveBlending,
            depthWrite: false,
            transparent: true,
        }));
        nearStarsGroup.add(sprite);
    }
    sprite.visible = true;
    sprite.material.opacity = opacity;
    sprite.material.color.set(getStarColor(star.spectralType));
    sprite.position.copy(starPosition);

    // フレアの大きさを計算 (距離と等級に依存)
    const scale = (starRenderParams.flareDistanceFactor / distance) * magnitudeFactor;
    const finalScale = THREE.MathUtils.clamp(scale, starRenderParams.flareMinScale, starRenderParams.flareMaxScale);
    sprite.scale.set(finalScale, finalScale, 1);
}


// --- 星のLOD（詳細度）を更新する ---
function updateStarLOD() {
    // 必要なデータがロードされるまで何もしない
    if (!starDataCache || !starParticles || !flareTexture || !originalSizes) return;

    const cameraPosition = camera.position;
    const pointSizes = starParticles.geometry.getAttribute('size');
    const positions = starParticles.geometry.getAttribute('position');
    let activeSphereCount = 0;
    let activeSpriteCount = 0;
    let needsSizeUpdate = false;
    const starPosition = new THREE.Vector3(); // Vector3オブジェクトを再利用して効率化

    // 全ての星についてカメラとの距離をチェック
    for (let i = 0; i < starDataCache.length; i++) {
        starPosition.fromBufferAttribute(positions, i);
        const distance = cameraPosition.distanceTo(starPosition);
        const star = starDataCache[i];
        const magnitude = star.magnitude;

        // 星の等級（明るさ）に基づく倍率を計算 (スプライトと球体で共通)
        const magnitudeFactor = Math.pow(
            starRenderParams.brightnessLogBase,
            (starRenderParams.brightnessReferenceMag - magnitude) / 2.0
        );

        // --- デフォルト状態を定義 ---
        let sphereOpacity = 0.0;
        let spriteOpacity = 0.0;
        let showParticle = true;

        // --- 距離に基づいて、表示するLODと透明度を決定 ---
        if (distance < lodDistance) {
            // 中景または最近景の範囲に入ったら、遠景のパーティクルは表示しない
            showParticle = false;

            // 遷移ゾーンの上限距離
            const transitionEnd = lodSphereDistance + lodTransitionRange;

            if (distance < lodSphereDistance) {
                // 最近景ゾーン: 球体のみ表示
                sphereOpacity = 1.0;
                spriteOpacity = 0.0;
            } else if (distance < transitionEnd) {
                // 遷移ゾーン: 球体とフレアをクロスフェード
                // progressが0(近い) -> 1(遠い)になるように計算
                const progress = (distance - lodSphereDistance) / lodTransitionRange;
                sphereOpacity = 1.0 - progress;
                spriteOpacity = progress;
            } else { // lodSphereDistance + lodTransitionRange <= distance < lodDistance
                // 中景ゾーン: フレアのみ表示
                sphereOpacity = 0.0;
                spriteOpacity = 1.0;
            }
        }

        // --- 決定した状態をオブジェクトに適用 ---

        // 1. パーティクルの表示/非表示を更新
        if (showParticle) {
            if (pointSizes.getX(i) === 0.0) {
                pointSizes.setX(i, originalSizes[i]);
                needsSizeUpdate = true;
            }
        } else {
            if (pointSizes.getX(i) > 0.0) {
                pointSizes.setX(i, 0.0);
                needsSizeUpdate = true;
            }
        }

        // 2. 球体の表示/透明度を更新
        if (sphereOpacity > 0) {
            updateLODSphere(activeSphereCount, star, starPosition, sphereOpacity, magnitudeFactor);
            activeSphereCount++;
        }

        // 3. スプライトの表示/透明度を更新
        if (spriteOpacity > 0) {
            updateLODSprite(activeSpriteCount, star, starPosition, spriteOpacity, magnitudeFactor, distance);
            activeSpriteCount++;
        }
    }

    // 使われなかったプール内の球体を非表示にする
    for (let i = activeSphereCount; i < nearSphereGroup.children.length; i++) {
        nearSphereGroup.children[i].visible = false;
    }

    // 使われなかったプール内のスプライトを非表示にする
    for (let i = activeSpriteCount; i < nearStarsGroup.children.length; i++) {
        nearStarsGroup.children[i].visible = false;
    }

    // サイズ属性に変更があった場合のみGPUにデータを転送
    if (needsSizeUpdate) {
        pointSizes.needsUpdate = true;
    }
}

// --- CSVテキストを解析する ---
function parseCSV(text) {
    const lines = text.trim().split('\n');
    const stars = [];
    let header = [];
    let headerParsed = false;
    const columnMap = {};

    for (const line of lines) {
        // コメント行と空行をスキップ
        if (line.startsWith('#') || line.trim() === '') {
            continue;
        }

        const tokens = line.split(',').map(t => t.trim());

        // 最初の非コメント行をヘッダーとして解析
        if (!headerParsed) {
            header = tokens;
            headerParsed = true;

            // ヘッダー名からインデックスを動的にマッピング
            columnMap.x = header.indexOf('x');
            columnMap.y = header.indexOf('y');
            columnMap.z = header.indexOf('z');
            columnMap.spectralType = header.indexOf('スペクトル');
            columnMap.magnitude = header.indexOf('星系視等級');

            // 必要な列が見つかるかチェック
            if (Object.values(columnMap).some(index => index === -1)) {
                console.error("CSVヘッダーエラー: 必要な列が見つかりません。", columnMap);
                throw new Error("CSV header is missing required columns (x, y, z, スペクトル, 星系視等級).");
            }
            continue; // ヘッダー行の処理はここまで
        }

        if (tokens.length < header.length) continue;

        const magnitude = parseFloat(tokens[columnMap.magnitude]);
        // 伴星（等級-1）や無効なデータは描画しない
        if (magnitude === -1.0 || isNaN(magnitude)) {
            continue;
        }

        const x = parseFloat(tokens[columnMap.x]);
        const y = parseFloat(tokens[columnMap.y]);
        const z = parseFloat(tokens[columnMap.z]);

        // 座標が無効な場合もスキップ
        if (isNaN(x) || isNaN(y) || isNaN(z)) continue;

        stars.push({
            x: x, y: y, z: z,
            spectralType: tokens[columnMap.spectralType],
            magnitude: magnitude
        });
    }
    return stars;
}

// --- 星のパーティクルシステムを生成する ---
async function createStars(starsData) {
    const geometry = new THREE.BufferGeometry();
    const positions = [];
    const colors = [];
    const sizes = [];
    const magnitudes = []; // 等級をシェーダーに渡すための配列

    const color = new THREE.Color();

    // ヘルパー関数: 等級からサイズを取得
    const getStarSize = (magnitude) => {
        // 等級が低い（明るい）ほど大きく表示。パラメータはstarRenderParamsで管理
        const size = starRenderParams.sizeMultiplier * Math.pow(starRenderParams.sizeLogBase, starRenderParams.sizeLogFactor * magnitude);
        return Math.max(starRenderParams.minSize, Math.min(size, starRenderParams.maxSize)); // サイズを適切な範囲にクランプ
    };

    for (const star of starsData) {
        positions.push(star.x, star.y, star.z);
        color.set(getStarColor(star.spectralType));
        colors.push(color.r, color.g, color.b);
        sizes.push(getStarSize(star.magnitude));
        magnitudes.push(star.magnitude); // 等級を追加
    }
    // LOD切り替えのために、元のサイズを保持しておく
    originalSizes = new Float32Array(sizes);

    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    // size属性は動的に変更するため、DynamicDrawUsageを設定
    const sizeAttribute = new THREE.Float32BufferAttribute(sizes, 1);
    sizeAttribute.setUsage(THREE.DynamicDrawUsage);
    geometry.setAttribute('size', sizeAttribute);
    geometry.setAttribute('magnitude', new THREE.Float32BufferAttribute(magnitudes, 1)); // 等級をAttributeとして設定

    // 必要なテクスチャを並行して非同期で読み込む
    const textureLoader = new THREE.TextureLoader();
    let pointTex;
    try {
        [pointTex, flareTexture] = await Promise.all([
            textureLoader.loadAsync('https://threejs.org/examples/textures/sprites/disc.png'),
            textureLoader.loadAsync('flare.png')
        ]);
    } catch (error) {
        console.error("Failed to load textures:", error);
        const info = document.getElementById('info');
        info.textContent = "Error: Could not load texture files. Check console for details.";
        info.style.color = "red";
        return; // テクスチャがないと進めないので処理を中断
    }

    // カスタムシェーダーを使ったマテリアル
    const material = new THREE.ShaderMaterial({
        uniforms: {
            pointTexture: { value: pointTex },
            uMinBrightness: { value: starRenderParams.minBrightness },
            uMaxBrightness: { value: starRenderParams.maxBrightness },
            uBrightnessLogBase: { value: starRenderParams.brightnessLogBase },
            uBrightnessReferenceMag: { value: starRenderParams.brightnessReferenceMag },
            uDistanceFalloffScale: { value: starRenderParams.distanceFalloffScale }
        },
        vertexShader: `
            attribute float size;
            attribute float magnitude;
            varying vec3 vColor;
            varying float vMagnitude;
            varying float vDistance;

            void main() {
                vColor = color;
                vMagnitude = magnitude;
                vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);

                vDistance = length(mvPosition.xyz);

                gl_PointSize = size * (300.0 / -mvPosition.z);
                gl_Position = projectionMatrix * mvPosition;
            }
        `,
        fragmentShader: `
            uniform sampler2D pointTexture;
            uniform float uMinBrightness;
            uniform float uMaxBrightness;
            uniform float uBrightnessLogBase;
            uniform float uBrightnessReferenceMag;
            uniform float uDistanceFalloffScale;
            varying vec3 vColor;
            varying float vMagnitude;
            varying float vDistance;

            void main() {
                vec4 discTex = texture2D(pointTexture, gl_PointCoord);

                if (discTex.a < 0.01) {
                    discard;
                }

                // 物理的な明るさの係数を計算 (等級が低いほど明るい)
                float brightnessMultiplier = pow(uBrightnessLogBase, uBrightnessReferenceMag - vMagnitude);
                float finalBrightness = clamp(brightnessMultiplier, uMinBrightness, uMaxBrightness);

                // カメラからの距離の2乗に応じて明るさを減衰させる
                // 近くの星が明るくなりすぎないように、減衰係数は最大1.0fにクランプする
                float distanceFalloff = clamp(uDistanceFalloffScale / (vDistance * vDistance), 0.0, 1.0);

                vec3 finalColor = vColor * finalBrightness * discTex.rgb * distanceFalloff;

                // AdditiveBlendingなのでアルファは1.0で良い
                gl_FragColor = vec4(finalColor, 1.0); // discTex.rgbでマスク済み
            }
        `,
        blending: THREE.AdditiveBlending,
        depthWrite: false, // 深度バッファへの書き込みはしない
        depthTest: true,   // 深度テストは行う
        transparent: true,
        vertexColors: true
    });

    starParticles = new THREE.Points(geometry, material);
    scene.add(starParticles);

    const info = document.getElementById('info');
    info.textContent = `${starsData.length} stars loaded. Use mouse to navigate.`;
}

// --- メイン処理の開始 ---
init();
animate();