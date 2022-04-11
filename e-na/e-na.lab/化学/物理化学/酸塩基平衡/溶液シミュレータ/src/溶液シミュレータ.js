class Util {
  static sum(...data) {
    return (Array.isArray(data[0]) ? data[0] : data).reduce((a, b) => a + b, 0);
  }

  static prod(...data) {
    return (Array.isArray(data[0]) ? data[0] : data).reduce((a, b) => a * b, 1);
  }
  static x(pX) {
    return Math.pow(10, -pX);
  }
  static p(x) {
    return Math.log10(x) * -1;
  }
    
  static flat(arr, d = 1) {
   return d > 0 ? arr.reduce((acc, val) => acc.concat(Array.isArray(val) ? flatDeep(val, d - 1) : val), [])
                : arr.slice();
  }
  static encodeHTML(formula, e){
    let words = formula.split(/[0-9]+/);
    let match = formula.match(/[0-9]+/g) || [];
      
    if(!words)return "";
      
    let array = [words[0]];

    for(let i = 0; i < match.length; ++ i){
      array.push(`<sub>${match[i]}</sub>`);
      array.push(words[i+1]);
    }

    array = array.join("");
    
    if(e != 0)array += `<sup>${(e==1||e==-1)?"":Math.abs(e)}${e>0?"+":"-"}</sup>`;
      
    console.log(array);
    return array;
  }
}

class Acid {
  //ほかのコードでこれだけ使うかも
  constructor(n = 1, pk = [0], charge = 0) {
    if (!Array.isArray(pk)) pk = [pk];

    this.n = n; //価数				int
    this.k = pk.map((pki) => Util.x(pki)); //酸解離定数			Array<double>[n]
    this.e = charge; //最多プロトン時の電荷		double
  }
}


class AcidInput extends Acid {
  //投入量付き
  constructor(n, pk, charge, c, name, formulas) {
    super(n, pk, charge);

    if (!Array.isArray(c)) c = [c]; //各イオンの投入量		Array<double>[n+1]
    this.c0 = c;

    this.c = this.c0.concat();
    this.name = name || `unknown${Date.now()}`;
    this.formulas = formulas;

    this.sum_c = Util.sum(this.c0);

    this.A = [1];

    for (let i = 1; i <= this.n; ++i)
      this.A.push(Util.prod(this.k.slice(0, i)));

    this.e0 = Util.sum(this.c.map((ci, i) => ci * (this.e - i)));

    //console.log(this);
  }

  setRatio(ph) {
    //各イオンの濃度[mol/L]を配列にして返す(i番目はi価(i = 0,1,...,n))
    let h = Util.x(ph);
    let sum = 0;

    for (let i = 0; i <= this.n; ++i) {
      this.c[i] = this.A[i] * Math.pow(h, this.n - i);
      sum += this.c[i];
    }

    this.c = this.c.map((ci) => (this.sum_c * ci) / sum);
    return this.c;
  }

  getCharge(ph) {
    //イオンによる1[L]当たりの電荷[C]をとってくる
    let h = Util.x(ph);
    let sum1 = 0,
      sum2 = 0;

    for (let i = 0; i <= this.n; ++i) {
      sum1 += this.A[i] * Math.pow(h, this.n - i);
      sum2 += this.A[i] * Math.pow(h, this.n - i) * (this.e - i);
    }

    //console.log(this.A);
    //console.log(sum1);

    return (this.sum_c * sum2) / sum1;
  }
}

class Acids {
  constructor(pksol = 14, data = []) {
    this.ksol = Util.x(pksol);

    if (!Array.isArray(data)) data = [data];
    this.acids = data.map((d) => new AcidInput(...d));

    this.sum_e0 = Util.sum(this.acids.map((acid) => acid.e0));

    //console.log(this);
  }

  test(ph) {
    let sum_e = Util.sum(this.acids.map((acid) => acid.getCharge(ph)));

    //console.log(sum_e + " at pH = " + ph);

    sum_e += Util.x(ph) - this.ksol / Util.x(ph);

    return sum_e - this.sum_e0;
  }

  get_pH() {
    let pH = 7;
    let d = 7 * (this.test(pH) < 0 ? -1 : 1);
    let max = (Math.log10(Math.abs(d)) + 16) / Math.log10(2);

    while (this.test(pH + d) * d > 0) {
      //console.log(`${this.test(pH + d)} at ${pH + d}`);
      pH += d;

      if (pH > 100) return -1;
    }

    for (let i = 0; i < max; ++i) {
      d /= 2;
      pH += d;

      d *= this.test(pH) * d < 0 ? -1 : 1;
    }
    
    pH += d;

    return [pH, Math.abs(d)];
  }

  setRatio() {
    let pHs = this.get_pH();

    this.pH = pHs[0];
    this.pH_diff = pHs[1];
    this.acids.forEach((acid) => acid.setRatio(pHs[0]));
  }

  property() {
    let str = "pKsol = " + Util.p(this.ksol) + "\n";
    let len = this.acids.length;
    let e = 0;
    let digit = Math.pow(10, 10);

    //str += `pH = ${this.pH}±${this.pH_diff}\n\n`;
    str += `pH = ${Math.round(this.pH*digit)/digit}\n\n`;

    this.acids.forEach((acid, i) => {
      console.log(acid);

      e = acid.e;
      str += `acid No.${i + 1}/${len}${acid.name?" ("+acid.name+")":""}:\n`;
      acid.c.forEach((ci, i) => {
        str += `  [${acid.formulas?Util.encodeHTML(acid.formulas[i], e-i):"charge="+(e-i)}]:${ci}←${acid.c0[i]}\n`;
      });

      str += "\n";
    });

    return str;
  }
}

class AcidPreset extends Acid {
  constructor(n, pk, charge, name, ion_formula) {
    super(n, pk, charge);
      
    this.name = name || `unknown${Date.now()}`;
    this.formulas = ion_formula;
  }
}

class AcidPresetCollecter {
  constructor(acidPresets){
    this.acidPresets = acidPresets;    
  }
    
  push(acidPreset){
    this.acidPresets.push(acidPreset); 
  }
}


function removeAcidPreset(target) {
  let index = target.parentNode.children[0].children[1].children[0];
    
  console.log(index);
  userPresets.acidPresets.splice(index, 1);

  target.parentNode.parentNode.remove();
}
function editAcidPreset(target) {
  let index = target.parentNode.children[0].children[0].innerHTML-1;
    
  console.log(index);
  userPresets.acidPresets.splice(index, 1, new AcidPreset(...get_preset_edit(target.parentNode.parentNode)));

  reloadPreset(target.parentNode.parentNode.parentNode.parentNode.parentNode.children[1]);
}
function registerPreset(target){
    userPresets.push(new AcidPreset(...get_preset_input(target.parentNode.children[0])));
    reloadPreset(target.parentNode.parentNode.parentNode.children[1]);
}
function resetPreset(target){
    userPresets = new AcidPresetCollecter([]);
    reloadPreset(target.parentNode.parentNode.parentNode.children[1]);
}
function get_preset_input(target){
    let n = 0, name = "", e = 0, pk = [], formula = [];
    
    name = target.children[0].children[0].children[0].value;
    e = target.children[1].children[0].children[0].value;
    formula.push(target.children[1].children[0].children[1].value);
    n = target.children[1].children[1].children.length;
    Array.prototype.forEach.call(target.children[1].children[1].children, ion =>{
        pk.push(ion.children[0].children[1].value);
        formula.push(ion.children[2].children[1].value);
    });
    
    console.log([n, pk, e, name, formula]);
    
    return [n, pk, e, name, formula];
}
function get_preset_edit(target){
    let n = 0, name = "", e = 0, pk = [], formula = [];
    
    console.log(target);
    
    name = target.children[0].children[0].children[1].value;
    e = target.children[1].children[0].children[0].value;
    formula.push(target.children[1].children[0].children[1].value);
    n = target.children[1].children[1].children.length;
    Array.prototype.forEach.call(target.children[1].children[1].children, ion =>{
        pk.push(ion.children[0].children[1].value);
        formula.push(ion.children[2].children[1].value);
    });
    
    console.log([n, pk, e, name, formula]);
    
    return [n, pk, e, name, formula];
}

function addPresetpKa(target) {
  let n = target.parentNode.children[1].childElementCount;

  target.parentNode.children[1].insertAdjacentHTML("beforeend", `<div class="ion">
                      <span class="pKa-in">
                        ↓pK<sub>a<span class="pka-sub">${n + 1}</span></sub> =
                        <input class="pKa" type="text" value="0" />
                      </span>

                      <br />

                      <span class="chemical">
                        電荷 <span class="charge">-1</span>:
                        <input class="formula" type="text" value="HSO4" />
                      </span>
                    </div>`);

  charge_reload(target.parentNode.children[0].children[0]);
}

function reloadPreset(target){
    target.parentNode.children[2].children[0].innerHTML = "";
    reloadDefaultPreset(target);
    reloadUserPreset(target);
    
    reloadSolution(target.parentNode.parentNode.parentNode.children[1].children[0]);
}
function reloadUserPreset(target){
    let table = target.parentNode.children[2].children[0];
    let html = "";
    
    userPresets.acidPresets.forEach((preset, index) => {
        html += `<div class="acid-preset-box">
                <span class="acid-preset-tab">
                  <h3 class="acid-preset-title">
                    User[<span id="index">${index+1}</span>]: <input type="text" class="preset-name" value="${preset.name}">
                  </h3>
                  <input type="button" class="remove-preset" value="削除" onclick="removeAcidPreset(this)">
                  <input type="button" class="edit-preset" value="適用" onclick="editAcidPreset(this)">
                </span>

                <div class="acid-preset-data">
                  <span class="chemical">
                    電荷<input type="text" class="preset-name" value="${preset.e}">: <input type="text" class="preset-name" value="${preset.formulas[0]}">(display: ${Util.encodeHTML(preset.formulas[0], preset.e)})
                  </span>
                  <div class="ions">
            `;
        
        for(let i = 0; i < preset.n; ++ i){
            html += `<div class="ion">
                      <span class="pKa-preset-in">
                        ↓pK<sub>a<span class="pka-preset-sub">${i+1}</span></sub> =
                        <input type="text" class="pka-preset" value="${Util.p(preset.k[i])}">
                      </span>

                      <br />

                      <span class="chemical">
                        電荷 <span class="charge">${preset.e-i-1}</span>:
                        <input type="text" class="preset-name" value="${preset.formulas[i+1]}">(display: ${Util.encodeHTML(preset.formulas[i+1], preset.e-i-1)})
                      </span>
                    </div>`;
        }
        
        html += `</div>
                 <input
                    type="button"
                    class="pKa-add"
                    value="電離定数の追加"
                    onclick="addPresetpKa(this);"
                  />
                  <input
                    type="button"
                    class="pKa-del"
                    value="電離定数の削除"
                    onclick="remove_pKa(this);"
                  />
                </div>
              </div>`;
    });
    
    //console.log(html);
    table.insertAdjacentHTML("beforeend", html);
}
function reloadDefaultPreset(target){
    let table = target.parentNode.children[2].children[0];
    let html = "";
    
    defaultPresets.acidPresets.forEach((preset, index) => {
        html += `<div class="acid-preset-box">
                <span class="acid-preset-tab">
                  <h3 class="acid-preset-title">
                    ${index+1}: ${preset.name}
                  </h3>
                </span>

                <div class="acid-preset-data">
                  <span class="chemical">
                    電荷${preset.e}: ${Util.encodeHTML(preset.formulas[0], preset.e)}
                  </span>
                  <div class="ions">
            `;
        
        for(let i = 0; i < preset.n; ++ i){
            html += `<div class="ion">
                      <span class="pKa-preset-in">
                        ↓pK<sub>a<span class="pka-preset-sub">${i+1}</span></sub> =
                        ${Util.p(preset.k[i])}
                      </span>

                      <br />

                      <span class="chemical">
                        電荷 <span class="charge">${preset.e-i-1}</span>:
                        ${Util.encodeHTML(preset.formulas[i+1], preset.e-i-1)}
                      </span>
                    </div>`;
        }
        
        html += `</div>
                </div>
              </div>`;
    });
    
    //console.log(html);
    table.insertAdjacentHTML("beforeend", html);
}

function get_sol_data(solution){
    return solution.children[1].children[1].value-0;
}
function get_acid_data(acids) {
  let acid_data = [],
    pKa = [],
    c0 = [];
  let id = "";

  let acids_data = [];

  //console.log(acids.children);

  Array.prototype.forEach.call(acids.children, (acid) => {
    acid_data = [acid.children[1].children[1].children.length];

    c0 = [acid.children[1].children[0].children[1].value - 0];

    pKa = Array.prototype.map.call(
      acid.children[1].children[1].children,
      (ion) => {
        c0.push(ion.children[2].children[1].value - 0);

        return ion.children[0].children[1].value - 0;
      }
    );

    acid_data.push(pKa);
    acid_data.push(acid.children[1].children[0].children[0].value - 0);

    acid_data.push(c0);
      
    id = acid.dataset.mode;
    if(id != -1){
        if(id.startsWith("d")){
            acid_data.push(defaultPresets.acidPresets[id.slice(1)].name);
            acid_data.push(defaultPresets.acidPresets[id.slice(1)].formulas);
        }
        else if(id.startsWith("u")){
            acid_data.push(userPresets.acidPresets[id.slice(1)].name);
            acid_data.push(userPresets.acidPresets[id.slice(1)].formulas);
        }
    }

    acids_data.push(acid_data);
  });

  console.log(acids_data);
  return acids_data;
}

function addAcid(target) {
  let acids = target.parentNode.children[0];
  let mode = target.parentNode.children[1].value;

  if(mode == -1)addAcidByHand(acids, mode);
  else if(mode.startsWith("d"))addAcidFromDefaultPreset(acids, mode.slice(1), mode);
  else if(mode.startsWith("u"))addAcidFromUserPreset(acids, mode.slice(1), mode);

  acids.lastChild.innerHTML += `<div class="preset-type" value=${mode}></div>`;
  acid_num_reload(acids);
}
function addAcidByHand(acids, mode) {
  acids.insertAdjacentHTML("beforeend", `<div class="acid-box" data-mode="${mode}">
                <span class="acid-tab">
                  <h3 class="acid-title">
                    酸 <span class="acid-num">1</span> of
                    <span class="acid-total">1</span>
                  </h3>
                  <input
                    type="button"
                    class="acid-del"
                    value="x"
                    onclick="removeAcid(this);"
                  />
                </span>

                <div class="acid-data">
                  <span class="chemical">
                    電荷
                    <input
                      class="charge"
                      type="text"
                      value="0"
                      onkeydown="charge_reload(this);"
                    />:
                    <input class="amount" type="text" value="0" />
                    M
                  </span>
                  <div class="ions"></div>
                  <input
                    type="button"
                    class="pKa-add"
                    value="電離定数の追加"
                    onclick="add_pKa(this);"
                  />
                  <input
                    type="button"
                    class="pKa-del"
                    value="電離定数の削除"
                    onclick="remove_pKa(this);"
                  />
                </div>
              </div>`);
}
function addAcidFromDefaultPreset(acids, i, mode) {
  let html = `<div class="acid-box" data-mode="${mode}">
                <span class="acid-tab">
                  <h3 class="acid-title">
                    酸 <span class="acid-num">1</span> of
                    <span class="acid-total">1</span>
                    (${defaultPresets.acidPresets[i].name})
                  </h3>
                  <input
                    type="button"
                    class="acid-del"
                    value="x"
                    onclick="removeAcid(this);"
                  />
                </span>

                <div class="acid-data">
                  <span class="chemical">
                    電荷
                    <input
                      class="charge"
                      type="text"
                      value="${defaultPresets.acidPresets[i].e}"
                      readonly
                    />:
                    <input class="amount" type="text" value="0" />
                    M     (${Util.encodeHTML(defaultPresets.acidPresets[i].formulas[0], defaultPresets.acidPresets[i].e)})
                  </span>
                  <div class="ions">`
  
         defaultPresets.acidPresets[i].k.forEach((kj, j) => {
             html += `<div class="ion">
                      <span class="pKa-in">
                        ↓pK<sub>a<span class="pka-sub">${j + 1}</span></sub> =
                        <input class="pKa" type="text" value="${Util.p(kj)}" readonly/>
                      </span>

                      <br />

                      <span class="chemical">
                        電荷 <span class="charge">${defaultPresets.acidPresets[i].e-j-1}</span>:
                        <input class="amount" type="text" value="0" />M     (${Util.encodeHTML(defaultPresets.acidPresets[i].formulas[j+1], defaultPresets.acidPresets[i].e-j-1)})
                      </span>
                    </div>\n`
         })
  
         html += `</div>
                </div>
              </div>`
  acids.insertAdjacentHTML("beforeend", html);
}
function addAcidFromUserPreset(acids, i, mode) {
    let html = `<div class="acid-box" data-mode="${mode}">
                <span class="acid-tab">
                  <h3 class="acid-title">
                    酸 <span class="acid-num">1</span> of
                    <span class="acid-total">1</span>
                    (${userPresets.acidPresets[i].name})
                  </h3>
                  <input
                    type="button"
                    class="acid-del"
                    value="x"
                    onclick="removeAcid(this);"
                  />
                </span>

                <div class="acid-data">
                  <span class="chemical">
                    電荷
                    <input
                      class="charge"
                      type="text"
                      value="${userPresets.acidPresets[i].e}"
                      readonly
                    />:
                    <input class="amount" type="text" value="0" />
                    M     (${Util.encodeHTML(userPresets.acidPresets[i].formulas[0], userPresets.acidPresets[i].e)})
                  </span>
                  <div class="ions">`
  
         userPresets.acidPresets[i].k.forEach((kj, j) => {
             html += `<div class="ion">
                      <span class="pKa-in">
                        ↓pK<sub>a<span class="pka-sub">${j + 1}</span></sub> =
                        <input class="pKa" type="text" value="${Util.p(kj)}" readonly/>
                      </span>

                      <br />

                      <span class="chemical">
                        電荷 <span class="charge">${userPresets.acidPresets[i].e-j-1}</span>:
                        <input class="amount" type="text" value="0" />M     (${Util.encodeHTML(userPresets.acidPresets[i].formulas[j+1], userPresets.acidPresets[i].e)})
                      </span>
                    </div>\n`
         })
  
         html += `</div>
                </div>
              </div>`
  acids.insertAdjacentHTML("beforeend", html);
}
function removeAcid(target) {
  let acids = target.parentNode.parentNode.parentNode;

  target.parentNode.parentNode.remove();
  acid_num_reload(acids);
}

function add_pKa(target) {
  let n = target.parentNode.children[1].childElementCount;

  target.parentNode.children[1].insertAdjacentHTML("beforeend", `<div class="ion">
                      <span class="pKa-in">
                        ↓pK<sub>a<span class="pka-sub">${n + 1}</span></sub> =
                        <input class="pKa" type="text" value="0" />
                      </span>

                      <br />

                      <span class="chemical">
                        電荷 <span class="charge">-1</span>:
                        <input class="amount" type="text" value="0" />M
                      </span>
                    </div>`);

  charge_reload(target.parentNode.children[0].children[0]);
}
function remove_pKa(target) {
  target.parentNode.children[1].children[
    target.parentNode.children[1].children.length - 1
  ].remove();
}

function addSolution(target) {
  let options = genOptions();
    
  console.log(options);
  target.parentNode.children[0].insertAdjacentHTML("beforeend", `<div class="solution">
          <span class="solution-tab">
            <h2 class="solution-title">
              溶液 <span class="solution-num">1</span> of
              <span class="solution-total">1</span>
            </h2>
            <input
              type="button"
              class="solution-del"
              value="x"
              onclick="removeSolution(this);"
            />
          </span>
          <span class="solution-data">
            pK<sub>sol</sub>=
            <input
              type="text"
              class="pKsol"
              value="14"
            />
          </span>
          <div class="in-box">
            <div class="acid-boxes"></div>

            <select name="selecter">
            ${options}
            </select>
            <input
              type="button"
              class="acid-add"
              value="酸の追加"
              onclick="addAcid(this);"
            />
            <input
              type="button"
              class="submit"
              value="適用"
              onclick="submit(this);"
            />
          </div>
          <div class="out-box">適用を押してね</div>
        </div>`);

  reloadSolution(target.parentNode.children[0]);
}
function removeSolution(target) {
  let sols = target.parentNode.parentNode.parentNode;
  target.parentNode.parentNode.remove();
  reloadSolution(sols);
}

function submit(target) {
  let sol = new Acids(get_sol_data(target.parentNode.parentNode), get_acid_data(target.parentNode.children[0]));
  let out = target.parentNode.parentNode.children[3];

  sol.setRatio();

  out.innerHTML = sol.property().replace(/\n/g, "<br>");
}

function acid_num_reload(target) {
  console.log(target);

  let n = target.childElementCount;

  Array.prototype.forEach.call(target.children, (acid, i) => {
    acid.children[0].children[0].children[0].textContent = i + 1;
    acid.children[0].children[0].children[1].textContent = n;
  });
}
function charge_reload(target) {
  let n = target.value;

  if (isNaN(n)) n = 0;
  Array.prototype.forEach.call(
    target.parentNode.parentNode.children[1].children,
    (ion, i) => {
      ion.children[2].children[0].textContent = n - i - 1;
    }
  );
}
function reloadSolution(target) {
  let n = target.childElementCount;

  Array.prototype.forEach.call(target.children, (sol, i) => {
    sol.children[0].children[0].children[0].textContent = i + 1;
    sol.children[0].children[0].children[1].textContent = n;
    sol.children[2].children[1].innerHTML = genOptions();
  });
}

function genOptions(){
  let options = "<option value=-1>手動で入力</option>";
  defaultPresets.acidPresets.forEach((preset, i) => {
      options += `<option value="d${i}">${preset.name}</option>\n`
  });
  userPresets.acidPresets.forEach((preset, i) => options += `<option value="u${i}">${preset.name}</option>\n`);
    
  return options;
}

function saveCookie(){
    document.cookie = "test=success";
}
function loadCookie(){
    console.log(document.cookie);
}



let defaultPresetSet = [];

defaultPresetSet.push(new AcidPreset(1, [-10], 0, "過塩素酸(水中)", ["HClO4", "ClO4"]));
defaultPresetSet.push(new AcidPreset(1, [-11], 0, "ヨウ化水素酸(水中)", ["HI", "I"]));
defaultPresetSet.push(new AcidPreset(1, [-9], 0, "臭化水素酸(水中)", ["HBr", "Br"]));
defaultPresetSet.push(new AcidPreset(1, [-7], 0, "塩酸(水中)", ["HCl", "Cl"]));
defaultPresetSet.push(new AcidPreset(2, [-2, 1.92], 0, "硫酸(水中)", ["H2SO4", "HSO4", "SO4"]));
defaultPresetSet.push(new AcidPreset(1, [-2], 0, "硝酸(水中)", ["HI", "I"]));
defaultPresetSet.push(new AcidPreset(2, [0, -14], 1, "水", ["H3O", "H2O", "OH"]));
defaultPresetSet.push(new AcidPreset(1, [1], 0, "塩素酸(水中)", ["HClO3", "ClO3"]));
defaultPresetSet.push(new AcidPreset(2, [1.81, 7.19], 0, "亜硫酸(水中)", ["H2SO3", "HSO3", "SO3"]));
defaultPresetSet.push(new AcidPreset(3, [2.12, 7.21, 12.67], 0, "リン酸(水中)", ["H3PO4", "H2PO4", "HPO4", "PO4"]));
defaultPresetSet.push(new AcidPreset(1, [3.45], 0, "フッ化水素酸(水中)", ["HF", "F"]));
defaultPresetSet.push(new AcidPreset(1, [3.75], 0, "蟻酸(水中)", ["HCOOH", "HCOO"]));
defaultPresetSet.push(new AcidPreset(1, [4.76], 0, "酢酸(水中)", ["CH3COOH", "CH3COO"]));
defaultPresetSet.push(new AcidPreset(1, [5.25], 1, "ピリジン(水中)", ["pyH", "py"]));
defaultPresetSet.push(new AcidPreset(2, [6.37, 10.32], 0, "炭酸(水中)", ["H2CO3", "HCO3", "CO3"]));
defaultPresetSet.push(new AcidPreset(2, [7.04, 19], 0, "硫化水素(水中)", ["H2S", "HS", "S"]));
defaultPresetSet.push(new AcidPreset(1, [9.14], 0, "ホウ酸(水中)", ["B(OH)3", "B(OH)4"]));
defaultPresetSet.push(new AcidPreset(1, [9.25], 1, "アンモニア(水中)", ["NH4", "NH3"]));
defaultPresetSet.push(new AcidPreset(1, [10.32], 0, "シアン化水素酸(水中)", ["HCN", "CN"]));


let defaultPresets = new AcidPresetCollecter(defaultPresetSet);
let userPresets = new AcidPresetCollecter([]);