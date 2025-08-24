// import {
//   WebGLRenderer,
//   Vector3,
//   Scene,
//   Quaternion,
//   CylinderGeometry,
//   AmbientLight,
//   InstancedMesh,
//   Object3D,
//   PointLight,
//   PerspectiveCamera,
//   Color,
//   SphereGeometry,
//   MeshLambertMaterial
// } from "https://cdn.skypack.dev/three@0.137";

// import { OrbitControls } from "https://cdn.skypack.dev/three-stdlib@2.8.5/controls/OrbitControls";

// // // === UI ELEMENTS ===
// // const addBtn = document.getElementById("addBtn");
// // const searchBar = document.getElementById("searchBar");
// // const sceneContainer = document.getElementById("sceneContainer");

// // addBtn.addEventListener("click", () => {
// //   searchBar.style.display =
// //     searchBar.style.display === "none" ? "block" : "none";
// //   searchBar.focus();
// // });

// // searchBar.addEventListener("keydown", (e) => {
// //   if (e.key === "Enter" && searchBar.value.toLowerCase() === "add") {
// //     createMoleculeScene();
// //     searchBar.value = "";
// //     searchBar.style.display = "none";
// //   }
// // });

// // === UI ELEMENTS ===
// const addBtn = document.getElementById("addBtn");
// const overlay = document.getElementById("overlay");
// const searchBar = document.getElementById("searchBar");
// const options = document.querySelectorAll("#options li");
// const sceneContainer = document.getElementById("sceneContainer");
// const overlayContent = document.querySelector(".overlay-content");
// // Show overlay when plus button is clicked
// addBtn.addEventListener("click", () => {
//   overlay.style.display = "flex";
//   searchBar.value = "";
//   searchBar.focus();
// });

// // Hide overlay when clicking outside content
// overlay.addEventListener("click", (e) => {
//   if (e.target === overlay) {
//     overlay.style.display = "none";
//   }
// });

// // Filter options as user types
// searchBar.addEventListener("input", () => {
//   const query = searchBar.value.toLowerCase();
//   options.forEach((opt) => {
//     opt.style.display = opt.textContent.toLowerCase().includes(query)
//       ? "block"
//       : "none";
//   });
// });

// // Allow Enter key for quick "add"
// searchBar.addEventListener("keydown", (e) => {
//   if (e.key === "Enter" && searchBar.value.toLowerCase() === "add") {
//     addSceneWrapper("custom");
//    // overlay.style.display = "none";
//   }
// });


// options.forEach(option => {
//   option.addEventListener("click", () => {
//     // trigger animation
//     overlayContent.classList.add("pop");

//     // remove class after animation so it can be reused
//     overlayContent.addEventListener("animationend", () => {
//       overlayContent.classList.remove("pop");
//     }, { once: true });
//   });
// });

// // Handle option clicks → add new scene
// options.forEach((opt) => {
//   opt.addEventListener("click", () => {
//     addSceneWrapper(opt.dataset.scene); 
//    // overlay.style.display = "none";
//   });
// });

// // === Wrapper that calls your existing function and adds remove button ===
// function addSceneWrapper(type = "default") {
//   const wrapper = document.createElement("div");
//   wrapper.className = "sceneBox";

//   // Create a canvas container for your scene
//   const canvasContainer = document.createElement("div");
//   wrapper.appendChild(canvasContainer);

//   // Call your scene creation function
//   createMoleculeScene();

//   // Label
//   const label = document.createElement("p");
//   label.textContent = type;
//   wrapper.appendChild(label);

//   // Remove button
//   // const removeBtn = document.createElement("button");
//   // removeBtn.className = "removeBtn";
//   // removeBtn.textContent = "Remove";
//   // removeBtn.addEventListener("click", (e) => {
//   //   console.log("Clicked remove for:", wrapper);  // 🟢 log the wrapper
//   //   console.log("Parent of removeBtn:", removeBtn.parentNode);
//   //   wrapper.remove();  // should remove the ENTIRE sceneBox
//   // });
//   // wrapper.appendChild(removeBtn);

//   // // Finally attach to sceneContainer
//   // //   sceneBox.appendChild(wrapper); good for keeping the search box

//   // sceneContainer.appendChild(wrapper);
// }





// // === FUNCTION TO CREATE A NEW MOLECULE INSTANCE ===
// function createMoleculeScene() {
//     // Wrapper div
//     // const box = document.createElement("div");
//     // box.className = "sceneBox";
//     // sceneContainer.appendChild(box);
   
//   // make a scene box wrapper
//   // const sceneBox = document.createElement("div");
//   // sceneBox.classList.add("sceneBox");

//   // // make a canvas
//   // const canvas = document.createElement("canvas");
//   // sceneBox.appendChild(canvas);

//   // // make a remove button
//   // const removeBtn = document.createElement("button");
//   // removeBtn.classList.add("removeBtn");
//   // removeBtn.textContent = "Delete";

//   // // when clicked, remove this whole sceneBox
//   // removeBtn.addEventListener("click", () => {
//   //   sceneBox.remove();
//   // });

//   // // add button *inside* sceneBox
//   // sceneBox.appendChild(removeBtn);

//   // // finally, add sceneBox to container
//   // document.getElementById("sceneContainer").appendChild(sceneBox);


//   const sceneBox = document.createElement("div");
//   sceneBox.classList.add("sceneBox");

//   // canvas first (scene on top)
//   //const canvas = document.createElement("canvas");
//   //sceneBox.appendChild(canvas);

//   // delete button AFTER canvas (below)
//   const removeBtn = document.createElement("button");
//   removeBtn.classList.add("removeBtn");
//   removeBtn.textContent = "Delete";
//   removeBtn.addEventListener("click", () => {
//     sceneBox.remove();
//   });
//   sceneBox.appendChild(removeBtn);

//   document.getElementById("sceneContainer").appendChild(sceneBox);


//     // Scene setup
//     const scene = new Scene();
//     scene.background = new Color("rgb(70, 161, 126)");

//     // Carbons
//     const carbonGeometry = new SphereGeometry(3, 32, 32);
//     const carbonMaterial = new MeshLambertMaterial({
//         color: new Color("rgb(107, 182, 87)")
//     });
//     //One instanced mesh for all carbons
//     const carbonPositions = [
//         [0, 0, 0],
//         [10, 0, 0],
//         [-10, 0, 0],
//         [0, 10, 0],
//         [0, -10, 0],
//         [0, 0, 10],
//     ];

//     const carbonMesh = new InstancedMesh(
//         carbonGeometry,
//         carbonMaterial,
//         carbonPositions.length
//     );

//     const dummy = new Object3D();
//         carbonPositions.forEach((pos, i) => {
//         dummy.position.set(...pos);
//         dummy.updateMatrix();
//         carbonMesh.setMatrixAt(i, dummy.matrix);
//     });
//     scene.add(carbonMesh);

//     // Bonds
//     const bondGeometry = new CylinderGeometry(0.5, 0.5, 1, 6, 1, true);
//     const bondMaterial = new MeshLambertMaterial({
//         color: new Color("rgb(87, 160, 182)")
//     });

//     const bondData = [];
//     for (let i = 0; i < carbonPositions.length; i++) {
//         for (let j = i + 1; j < carbonPositions.length; j++) {
//             const start = new Vector3(...carbonPositions[i]);
//             const end = new Vector3(...carbonPositions[j]);
//             const midPoint = new Vector3().addVectors(start, end).multiplyScalar(0.5);
//             const direction = new Vector3().subVectors(end, start);
//             const length = direction.length();

//             const quaternion = new Quaternion();
//             quaternion.setFromUnitVectors(
//             new Vector3(0, 1, 0),
//             direction.clone().normalize()
//             );

//             bondData.push({ position: midPoint.toArray(), quaternion, length });
//         }
//     }
//     //One instanced mesh for all carbons
//     const bondMesh = new InstancedMesh(bondGeometry, bondMaterial, bondData.length);
//     const dummyBond = new Object3D();
//     bondData.forEach((bond, i) => {
//         dummyBond.position.set(...bond.position);
//         dummyBond.quaternion.copy(bond.quaternion);// orientation
//         dummyBond.scale.set(1, bond.length, 1);// stretch cylinder along Y-axis
//         dummyBond.updateMatrix();
//         bondMesh.setMatrixAt(i, dummyBond.matrix);
//     });
//     scene.add(bondMesh);

//     const ambientLight = new AmbientLight(new Color("rgb(255, 255, 255)"),1);//color, intensity
//     const pointLight = new PointLight(new Color("rgb(255, 255, 255)"),1,0,2);//color, intensity,distance (default 0 infinite),decay (default 2)
//     pointLight.position.set(20,20,0);
//     scene.add(ambientLight);
//     scene.add(pointLight);

//     const width = 250;
//     const height = 350;
//     const aspect = width / height; 
//     const camera = new PerspectiveCamera(45,aspect,0.1,1000);//check variables and decide ortho?  original (45,innerWidth/innerHeight,0.1,1000)
//     camera.position.set(0,0,50/aspect);

//     //remeber to change it in a way that it makes it compatible for mobile phones
//     const renderer = new WebGLRenderer({ antialias: true });
//     renderer.setPixelRatio(window.devicePixelRatio); //can divide by some number greater than 1 to boost performance
//     //renderer.setSize(window.innerWidth, window.innerHeight);
//     //renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
//     renderer.setSize(width, height);
//     sceneBox.appendChild(renderer.domElement);


//     const controls = new OrbitControls(camera,renderer.domElement);
//     controls.target.set(0,0,0); //always looks at origin
//     controls.enableDamping=true;//start is not dampened but end is
//     controls.dampingFactor= 0.03; //dampness of controls check later
//     controls.minPolarAngle = Math.PI / 2-0.25; // lock at 90°
//     controls.maxPolarAngle = Math.PI / 2+0.25; // lock at 90°
//     controls.enableZoom = false;
//     controls.enablePan = false;

//     renderer.setAnimationLoop(()=>{

//     controls.update();//need to call this every frame if damping or auto rotate is used
//     renderer.render(scene,camera);
    
//     });
// }



import {
  WebGLRenderer,
  Vector3,
  Scene,
  Quaternion,
  CylinderGeometry,
  AmbientLight,
  InstancedMesh,
  Object3D,
  PointLight,
  PerspectiveCamera,
  Color,
  SphereGeometry,
  MeshLambertMaterial,
  Mesh
} from "https://cdn.skypack.dev/three@0.137";

import { OrbitControls } from "https://cdn.skypack.dev/three-stdlib@2.8.5/controls/OrbitControls";

/* =========================
   RDKit loader (optional)
   ========================= */
let RDKit = null;
let rdkitReady = (async () => {
  try {
    // If RDKit script tag is present, init it
    if (typeof initRDKitModule === "function") {
      RDKit = await initRDKitModule();
      console.log("[RDKit] loaded");
    } else {
      console.warn("[RDKit] script tag not found. Falling back to 2D/placeholder coords.");
    }
  } catch (e) {
    console.warn("[RDKit] failed to init:", e);
  }
})();

/* =========================
   UI ELEMENTS
   ========================= */
const addBtn = document.getElementById("addBtn");
const overlay = document.getElementById("overlay");
const searchBar = document.getElementById("searchBar");
const options = document.querySelectorAll("#options li");
const sceneContainer = document.getElementById("sceneContainer");
const overlayContent = document.querySelector(".overlay-content");

// Show overlay when plus button is clicked
addBtn.addEventListener("click", () => {
  overlay.style.display = "flex";
  searchBar.value = "";
  searchBar.focus();
});

// Hide overlay when clicking outside content
overlay.addEventListener("click", (e) => {
  if (e.target === overlay) overlay.style.display = "none";
});

// Filter options as user types
searchBar.addEventListener("input", () => {
  const query = searchBar.value.toLowerCase();
  options.forEach((opt) => {
    opt.style.display = opt.textContent.toLowerCase().includes(query)
      ? "block"
      : "none";
  });
});

// Enter in searchBar now treats value as a SMILES string
searchBar.addEventListener("keydown", async (e) => {
  if (e.key === "Enter") {
    const smiles = searchBar.value.trim();
    if (!smiles) return;
    await addSceneWrapper(smiles);
    // overlay.style.display = "none"; // keep open if you want to add many quickly
    searchBar.select();
  }
});

// Click on list options → use the label as SMILES for quick testing
options.forEach(option => {
  option.addEventListener("click", async () => {
    overlayContent.classList.add("pop");
    overlayContent.addEventListener("animationend", () => {
      overlayContent.classList.remove("pop");
    }, { once: true });

    const smiles = option.dataset.scene ; // example fallback
    await addSceneWrapper(smiles);
    
  });
});

/* =========================
   Add Scene Wrapper
   ========================= */
async function addSceneWrapper(smiles = "CCO") {
  const sceneBox = document.createElement("div");
  sceneBox.className = "sceneBox";

  // Delete button
  const removeBtn = document.createElement("button");
  removeBtn.className = "removeBtn";
  removeBtn.textContent = "Delete";
  //removeBtn.addEventListener("click", () => sceneBox.remove());
removeBtn.addEventListener("click", () => {
  sceneBox.style.display = "none";   // hides it
  console.log(smiles +" deleted");
});
  sceneContainer.appendChild(sceneBox);

  // build the actual viewer
  await createMoleculeSceneFromSmiles(smiles, sceneBox);

  // place the button after canvas
  sceneBox.appendChild(removeBtn);
}

/* =========================
   SMILES → atoms/bonds
   =========================
   Returns:
   {
     atoms: [{ index, symbol, x, y, z }],
     bonds: [{ begin, end, order }]  // order: 1|2|3
   }
*/
async function smilesToMolData(smiles) {
  await rdkitReady;

  // Helper to parse a V2000/V3000 molblock (2D) into coords/bonds
  function parseMolblock(molblock) {
    const atoms = [];
    const bonds = [];

    if (!molblock) return { atoms, bonds };

    const lines = molblock.split(/\r?\n/);
    // counts line is line 4 (index 3) in V2000; V3000 is different.
    // We'll try V2000 path first and fallback to very basic V3000 parsing.
    const countsLine = lines[3] || "";
    const v2000Match = countsLine.match(/^\s*(\d+)\s+(\d+)/);
    let atomCount = 0, bondCount = 0, cursor = 4;

    if (v2000Match) {
      atomCount = parseInt(v2000Match[1]);
      bondCount = parseInt(v2000Match[2]);

      // atoms (x y z symbol ...)
      for (let i = 0; i < atomCount; i++) {
        const L = lines[cursor + i] || "";
        const x = parseFloat(L.slice(0, 10));
        const y = parseFloat(L.slice(10, 20));
        const z = parseFloat(L.slice(20, 30));
        const symbol = L.slice(31, 34).trim() || "C";
        atoms.push({ index: i, symbol, x, y, z: isFinite(z) ? z : 0 });
      }
      cursor += atomCount;

      // bonds (begin end order ...)
      for (let j = 0; j < bondCount; j++) {
        const L = lines[cursor + j] || "";
        const a1 = parseInt(L.slice(0, 3)) - 1;
        const a2 = parseInt(L.slice(3, 6)) - 1;
        const order = parseInt(L.slice(6, 9)) || 1;
        if (Number.isInteger(a1) && Number.isInteger(a2)) {
          bonds.push({ begin: a1, end: a2, order });
        }
      }
      return { atoms, bonds };
    }

    // extremely light V3000 fallback: look for ATOM/BOND blocks
    let inAtoms = false, inBonds = false;
    for (const L of lines) {
      const line = L.trim();
      if (line.includes("BEGIN ATOM")) { inAtoms = true; continue; }
      if (line.includes("END ATOM")) { inAtoms = false; continue; }
      if (line.includes("BEGIN BOND")) { inBonds = true; continue; }
      if (line.includes("END BOND")) { inBonds = false; continue; }

      if (inAtoms) {
        // e.g. "M  V30 1 C  -0.7500 0.0000 0.0000 0"
        const parts = line.split(/\s+/);
        const idx = parseInt(parts[2]) - 1;
        const symbol = parts[3] || "C";
        const x = parseFloat(parts[4]);
        const y = parseFloat(parts[5]);
        const z = parseFloat(parts[6]) || 0;
        atoms.push({ index: idx, symbol, x, y, z });
      } else if (inBonds) {
        // e.g. "M  V30 1 1  1 2"
        const parts = line.split(/\s+/);
        const order = parseInt(parts[3]) || 1;
        const a1 = parseInt(parts[4]) - 1;
        const a2 = parseInt(parts[5]) - 1;
        bonds.push({ begin: a1, end: a2, order });
      }
    }
    // normalize by index
    atoms.sort((a,b)=>a.index-b.index).forEach((a,i)=>a.index=i);
    return { atoms, bonds };
  }

  // If RDKit available, try to get 3D coords
  if (RDKit) {
    try {
      const mol = RDKit.get_mol(smiles);

      // Try to embed + optimize (3D). These names vary by RDKit.js build;
      // we attempt several and fall back to 2D molblock if needed.
      let have3D = false;

      // Attempt #1: built-in helper (if exposed)
      try {
        if (typeof mol.EmbedMolecule === "function") {
          mol.EmbedMolecule(); // ETKDG default
          have3D = true;
        }
      } catch {}

      // Attempt #2: force MMFF optimization if present
      try {
        if (typeof mol.MMFFOptimizeMolecule === "function") {
          mol.MMFFOptimizeMolecule();
        }
      } catch {}

      // Attempt to read coords (if API exposed)
      if (have3D && typeof mol.get_conformer === "function") {
        const atoms = [];
        const bonds = [];
        const nA = mol.get_num_atoms();
        const nB = mol.get_num_bonds();

        const conf = mol.get_conformer();
        for (let i = 0; i < nA; i++) {
          const atom = mol.get_atom_with_idx(i);
          const pos = conf.get_atom_position(i);
          atoms.push({
            index: i,
            symbol: atom.get_symbol(),
            x: pos.x, y: pos.y, z: pos.z
          });
        }
        for (let j = 0; j < nB; j++) {
          const b = mol.get_bond_with_idx(j);
          // get_bond_type() may be "SINGLE","DOUBLE","TRIPLE" in some builds
          let order = 1;
          try {
            const t = b.get_bond_type?.();
            if (t === "DOUBLE" || t === 2) order = 2;
            if (t === "TRIPLE" || t === 3) order = 3;
          } catch {}
          bonds.push({
            begin: b.get_begin_atom_idx(),
            end: b.get_end_atom_idx(),
            order
          });
        }
        mol.delete?.();
        if (atoms.length) return { atoms, bonds };
      }

      // Fallback: 2D molblock (z=0)
      const molblock = mol.get_molblock?.() || "";
      mol.delete?.();
      const data2D = parseMolblock(molblock);
      if (data2D.atoms.length) return data2D;
    } catch (e) {
      console.warn("[RDKit] SMILES→3D failed. Using fallback 2D:", e);
    }
  }

  // Final fallback: trivial linear layout for quick visualization
  const atoms = [];
  const bonds = [];
  const symbols = smiles.replace(/[^A-Z]/gi, "").match(/[A-Z][a-z]?/g) || ["C","C"];
  symbols.forEach((sym, i) => atoms.push({ index: i, symbol: sym, x: i*2.2, y: 0, z: 0 }));
  for (let i=0;i<symbols.length-1;i++) bonds.push({ begin: i, end: i+1, order: 1 });
  return { atoms, bonds };
}

/* =========================
   THREE: build scene
   ========================= */
async function createMoleculeSceneFromSmiles(smiles, mountEl) {
  const scene = new Scene();
  scene.background = new Color("rgb(70, 161, 126)");

  const width = 250;
  const height = 350;
  const aspect = width / height;

  const camera = new PerspectiveCamera(45, aspect, 0.1, 1000);
  camera.position.set(0, 0, 50 / aspect);

  const renderer = new WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(window.devicePixelRatio);
  renderer.setSize(width, height);
  mountEl.prepend(renderer.domElement); // canvas on top

  const controls = new OrbitControls(camera, renderer.domElement);
  controls.target.set(0, 0, 0);
  controls.enableDamping = true;
  controls.dampingFactor = 0.03;
  controls.minPolarAngle = Math.PI / 2 - 0.25;
  controls.maxPolarAngle = Math.PI / 2 + 0.25;
  controls.enableZoom = false;
  controls.enablePan = false;

  const ambientLight = new AmbientLight(new Color("rgb(255, 255, 255)"), 1);
  const pointLight = new PointLight(new Color("rgb(255, 255, 255)"), 1, 0, 2);
  pointLight.position.set(20, 20, 0);
  scene.add(ambientLight, pointLight);

  // Get atom/bond data
  const { atoms, bonds } = await smilesToMolData(smiles);

  // Center & scale the molecule nicely
  centerAtomPositions(atoms, 14);

  // --- ATOMS (spheres) ---
  const atomRadius = 1;
  const atomGeom = new SphereGeometry(atomRadius, 32, 32);

  // Instanced by element for color batching
  const byElement = groupBy(atoms, a => a.symbol);
  Object.entries(byElement).forEach(([symbol, list]) => {
    const mat = new MeshLambertMaterial({ color: cpkColor(symbol) });
    const mesh = new InstancedMesh(atomGeom, mat, list.length);
    const dummy = new Object3D();
    list.forEach((a, i) => {
      dummy.position.set(a.x, a.y, a.z);
      dummy.updateMatrix();
      mesh.setMatrixAt(i, dummy.matrix);
    });
    mesh.instanceMatrix.needsUpdate = true;
    scene.add(mesh);
  });

  // --- BONDS (cylinders) ---
  // We’ll draw one cylinder per *segment*.
  // Double/triple bonds are drawn as parallel cylinders offset slightly.
  const baseBondGeom = new CylinderGeometry(0.45, 0.45, 1, 10, 1, true);
  const bondMat = new MeshLambertMaterial({ color: new Color("rgb(34, 63, 72)") });

  const bondOffset = 0.6; // spacing for multi-bonds
  for (const b of bonds) {
    const a1 = atoms[b.begin];
    const a2 = atoms[b.end];
    if (!a1 || !a2) continue;

    const start = new Vector3(a1.x, a1.y, a1.z);
    const end   = new Vector3(a2.x, a2.y, a2.z);
    const dir = new Vector3().subVectors(end, start);
    const len = dir.length();
    if (len < 1e-6) continue;

    const mid = new Vector3().addVectors(start, end).multiplyScalar(0.5);
    const quat = new Quaternion().setFromUnitVectors(new Vector3(0, 1, 0), dir.clone().normalize());

    const segments = bondSegments(b.order, start, end, bondOffset);
    for (const seg of segments) {
      const segMid = new Vector3().addVectors(seg.start, seg.end).multiplyScalar(0.5);
      const segDir = new Vector3().subVectors(seg.end, seg.start);
      const segLen = segDir.length();
      const segQuat = new Quaternion().setFromUnitVectors(new Vector3(0, 1, 0), segDir.clone().normalize());

      const cyl = new Mesh(baseBondGeom, bondMat);
      cyl.position.copy(segMid);
      cyl.quaternion.copy(segQuat);
      cyl.scale.set(1, segLen, 1);
      scene.add(cyl);
    }
  }

  controls.autoRotate = true;      // turn on auto-rotation
controls.autoRotateSpeed = -10.0;  // default is 2.0, you can make it slower/faster

  renderer.setAnimationLoop(() => {
    controls.update();
    renderer.render(scene, camera);
  });
}

/* =========================
   Utils
   ========================= */

// center & scale molecule to fit a target radius
function centerAtomPositions(atoms, targetRadius = 12) {
  if (!atoms.length) return;
  const cx = atoms.reduce((s, a) => s + a.x, 0) / atoms.length;
  const cy = atoms.reduce((s, a) => s + a.y, 0) / atoms.length;
  const cz = atoms.reduce((s, a) => s + a.z, 0) / atoms.length;

  let maxR = 1e-6;
  atoms.forEach(a => {
    a.x -= cx; a.y -= cy; a.z -= cz;
    const r = Math.hypot(a.x, a.y, a.z);
    if (r > maxR) maxR = r;
  });

  const s = targetRadius / maxR;
  atoms.forEach(a => { a.x *= s; a.y *= s; a.z *= s; });
}

function groupBy(arr, keyFn) {
  return arr.reduce((m, item) => {
    const k = keyFn(item);
    (m[k] ||= []).push(item);
    return m;
  }, {});
}

// very small CPK palette
function cpkColor(symbol) {
  const table = {
    H: 0xFFFFFF,
    C: 0x222222,
    N: 0x3050F8,
    O: 0xFF0D0D,
    F: 0x90E050,
    Cl: 0x1FF01F,
    Br: 0xA62929,
    I: 0x940094,
    S: 0xFFFF30,
    P: 0xFF8000,
    B: 0xFFB5B5,
    Si: 0xBFBFBF
  };
  return table[symbol] ?? 0x888888;
}

// create parallel segments for multi-bonds
function bondSegments(order, start, end, offsetDist) {
  const segments = [];
  if (order <= 1) {
    segments.push({ start, end });
    return segments;
  }

  // direction and a perpendicular basis
  const dir = new Vector3().subVectors(end, start).normalize();
  // pick an arbitrary vector not parallel to dir
  const tmp = Math.abs(dir.y) < 0.9 ? new Vector3(0, 1, 0) : new Vector3(1, 0, 0);
  const perp = new Vector3().crossVectors(dir, tmp).normalize();
  const offset = perp.multiplyScalar(offsetDist);

  if (order === 2) {
    segments.push({ start: new Vector3().addVectors(start, offset), end: new Vector3().addVectors(end, offset) });
    segments.push({ start: new Vector3().subVectors(start, offset), end: new Vector3().subVectors(end, offset) });
  } else if (order >= 3) {
    // center + ±offset*1.2
    segments.push({ start, end });
    const off = offset.clone().multiplyScalar(1.2);
    segments.push({ start: new Vector3().addVectors(start, off), end: new Vector3().addVectors(end, off) });
    segments.push({ start: new Vector3().subVectors(start, off), end: new Vector3().subVectors(end, off) });
  }
  return segments;
}

/* =========================
   Drug Interaction Logic
   ========================= */

let selectedSmiles = [];

// 🔹 Handles selection of drug SMILES
function handleDrugSelection(smiles) {
  selectedSmiles.push(smiles);
  if (selectedSmiles.length === 2) {
    sendInteractionRequest(selectedSmiles[0], selectedSmiles[1]);
    selectedSmiles = []; // reset for next pair
  }
}

// 🔹 Sends POST request to Flask backend
async function sendInteractionRequest(smiles1, smiles2) {
  try {
    const response = await fetch("/predict", {
      method: "POST",
      headers: {
        "Content-Type": "application/json"
      },
      body: JSON.stringify({
        drug1_name: smiles1,
        drug2_name: smiles2
      })
    });

    const result = await response.json();
    if (response.ok) {
      showInteractionResult(result.prediction, smiles1, smiles2);
    } else {
      showInteractionError(result.error || "Prediction failed.");
    }
  } catch (err) {
    console.error("Interaction request failed:", err);
    showInteractionError("Something went wrong.");
  }
}

// 🔹 Displays result in a container
function showInteractionResult(prediction, smiles1, smiles2) {
  const container = document.getElementById("sceneContainer");
  container.innerHTML = `
    <div class="result">
      <h3>Interaction Result</h3>
      <p><strong>${smiles1}</strong> + <strong>${smiles2}</strong> → 
      ${prediction === 1 ? "⚠️ Interaction Detected" : "✅ No Interaction"}</p>
    </div>
  `;
}

// 🔹 Displays error messages
function showInteractionError(message) {
  const container = document.getElementById("sceneContainer");
  container.innerHTML = `
    <div class="error">
      <h3>Error</h3>
      <p>${message}</p>
    </div>
  `;
}

// 🔹 Attach click handlers to drug list items
document.addEventListener("DOMContentLoaded", () => {
  const drugList = document.getElementById("drugList");
  if (drugList) {
    drugList.addEventListener("click", (e) => {
      if (e.target && e.target.tagName === "LI") {
        const smiles = e.target.getAttribute("data-smiles");
        if (smiles) {
          handleDrugSelection(smiles);
        }
      }
    });
  }
});
