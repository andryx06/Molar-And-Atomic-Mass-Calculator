"use strict";

const atomicTable={
    "source": "http://www.chemistry.or.jp/activity/atomictable2020.pdf#page=4",
    "H":1.008,
    "He":4.003,
    "Li":6.941,
    "Be":9.012,
    "B":10.81,
    "C":12.01,
    "N":14.01,
    "O":16.00,
    "F":19.00,
    "Ne":20.18,
    "Na":22.99,
    "Mg":24.31,
    "Al":26.98,
    "Si":28.09,
    "P":30.97,
    "S":32.07,
    "Cl":35.45,
    "Ar":39.95,
    "K":39.10,
    "Ca":40.08,
    "Sc":44.96,
    "Ti":47.87,
    "V":50.94,
    "Cr":52.00,
    "Mn":54.94,
    "Fe":55.85,
    "Co":58.93,
    "Ni":58.69,
    "Cu":63.55,
    "Zn":65.38,
    "Ga":69.72,
    "Ge":72.63,
    "As":74.92,
    "Se":78.97,
    "Br":79.90,
    "Kr":83.80,
    "Rb":85.47,
    "Sr":87.62,
    "Y":88.91,
    "Zr":91.22,
    "Nb":92.91,
    "Mo":95.95,
    "Tc":99,
    "Ru":101.1,
    "Rh":102.9,
    "Pd":106.4,
    "Ag":107.9,
    "Cd":112.4,
    "In":114.8,
    "Sn":118.7,
    "Sb":121.8,
    "Te":127.6,
    "I":126.9,
    "Xe":131.3,
    "Cs":132.9,
    "Ba":137.3,
    "La":138.9,
    "Ce":140.1,
    "Pr":140.9,
    "Nd":144.2,
    "Pm":145,
    "Sm":150.4,
    "Eu":152.0,
    "Gd":157.3,
    "Tb":158.9,
    "Dy":162.5,
    "Ho":164.9,
    "Er":167.3,
    "Tm":168.9,
    "Yb":173.0,
    "Lu":175.0,
    "Hf":178.5,
    "Ta":180.9,
    "W":183.8,
    "Re":186.2,
    "Os":190.2,
    "Ir":192.2,
    "Pt":195.1,
    "Au":197.0,
    "Hg":200.6,
    "Tl":204.4,
    "Pb":207.2,
    "Bi":209.0,
    "Po":210,
    "At":210,
    "Rn":222,
    "Fr":223,
    "Ra":226,
    "Ac":227,
    "Th":232.0,
    "Pa":231.0,
    "U":238.0,
    "Np":237,
    "Pu":239,
    "Am":243,
    "Cm":247,
    "Bk":247,
    "Cf":252,
    "Es":252,
    "Fm":257,
    "Md":258,
    "No":259,
    "Lr":262,
    "Rf":267,
    "Db":268,
    "Sg":271,
    "Bh":272,
    "Hs":277,
    "Mt":276,
    "Ds":281,
    "Rg":280,
    "Cn":285,
    "Nh":278,
    "Fl":289,
    "Mc":289,
    "Lv":293,
    "Ts":293,
    "Og":294
};
const atomicWeightSource=document.getElementById('atomicWeight');
atomicWeightSource.href=atomicTable.source;

window.onload=()=>{calcCompound(); calcEquation();}

const braL=/[[({]/;
const braR=/[\])}]/;
const compoundFormula=document.getElementById('compoundFormula');
const outputCompoundWeight=document.getElementById('outputCompoundWeight');
const chemicalEquation=document.getElementById('chemicalEquation');
const outputChemicalEquation=document.getElementById('outputChemicalEquation');

// for compound
const calcCompound=()=>{
    outputCompoundWeight.value="";
    outputCompoundWeight.value=calcFormulaWeight("("+compoundFormula.value.replace(/\./,")(")+")");
};
const calcFormulaWeight=formula=>{
    const len=formula.length;
    let sumWeight=0;
    let iStart=0;
    let coef;
    if((coef=/^[0-9]+/.exec(formula))==null){
        coef=1;
    }else{
        iStart=coef.toString().length;
        coef=parseInt(coef);
    }
    for(let i=iStart;i<len;i++){
        let atom="";
        if(/[A-Z]/.test(formula[i])){
            atom+=formula[i];
            while(/[a-z0-9]/.test(formula[i+1])&&i+1<len) atom+=formula[1+i++];
            sumWeight+=decompose2Atoms(atom);
        }else if(braL.test(formula[i])){
            let nBra=1;
            while(nBra>0&&i+1<len){
                if(braR.test(formula[i+1])){
                    nBra--;
                    if(nBra>0) atom+=formula[i+1];
                }else if(braL.test(formula[i+1])){
                    nBra++;
                    atom+=formula[i+1];
                }else{
                    atom+=formula[i+1];
                }
                i++;
            }
            let coef2="";
            while(/[0-9]/.test(formula[i+1])&&i+1<len) coef2+=formula[1+i++];
            coef2=(coef2==="" ? 1 : parseInt(coef2));
            sumWeight+=coef2*calcFormulaWeight(atom);
        }else{
            console.log("undefined:",formula[i]);
        }
    }
    return coef*sumWeight;
};
const decompose2Atoms=(formula)=>{
    const atoms=formula.split(/([A-Z][a-z]?\d*)/).filter(e => e!=="");
    const coef=(isFinite(atoms[0]) ? atoms[0] : 1);
    return coef*atoms.slice(coef>1 ? 1 : 0).reduce((sum,e)=>sum+calcWeight(e.split(/(\d+)/).filter(e=>e!=="")),0);
};
const calcWeight=(atom)=>(atom.length===2?atomicTable[atom[0]]*atom[1]:atomicTable[atom]);

// for chemical equation
const calcEquation=()=>{
    outputChemicalEquation.value="";
    outputChemicalEquation.value=calcEquationCoef(chemicalEquation.value.replace(/\s*[+,]\s*/g,","));
};
const calcEquationCoef=(equation)=>{
    let ans="";
    let i=0,atom,tmp;
    let atoms=new Map();
    equation.split(/([A-Z][a-z]?[0-9]*)/).filter(e=>/[A-Za-z]+[0-9]*/.test(e)).forEach(e=>{
        atom=/[A-Z][a-z]?/.exec(e)[0];
        if(!atoms.has(atom)) atoms.set(atom,i++);
    });
    let [left,right]=equation.split(/=|=>|->|⇒|→/).map(e=>e.trim().split(","));
    let coef=new Array(atoms.size);
    for(i=0;i<coef.length;i++) coef[i]=new Array(right.length+left.length).fill(0);
    for(i=0;i<left.length;i++) setAtomCoef(i,"("+left[i].replace(/\./,")(")+")",atoms,coef,1);
    for(i=0;i<right.length;i++) setAtomCoef(left.length+i,"("+right[i].replace(/\./,")(")+")",atoms,coef,-1);
    pivotGaussJordan(coef);
    for(i=0;i<left.length;i++){
        ans+=((tmp=Math.round(coef[i][i]*1000)/1000)===1?"":tmp)+left[i]+" + ";
    }
    ans=ans.slice(0,-3)+" = ";
    for(i=left.length;i<left.length+right.length-1;i++){
        ans+=((tmp=-Math.round(coef[i][i]*1000)/1000)===1?"":tmp)+right[i-left.length]+" + ";
    }
    ans+=((tmp=Math.round(coef[0].pop()*1000)/1000)===1?"":tmp)+right.pop();
    return ans;
};

const setAtomCoef=(row,str,atoms,coefMatrix,coef)=>{
    let i,coef1,iStart=0;
    if((coef1=/^[0-9]+/.exec(str))==null){
        coef1=1;
    }else{
        iStart=coef1.toString().length;
        coef1=parseInt(coef1);
    }
    for(i=iStart;i<str.length;i++){
        let atom="";
        if(/[A-Z]/.test(str[i])){
            atom+=str[i];
            while(/[a-z0-9]/.test(str[i+1])&&i+1<str.length) atom+=str[1+i++];
            atom=atom.split(/([A-Z][a-z]?)/).filter(e=>/\S/.test(e));
            if(atom.length===1) atom.push("1");
            coefMatrix[atoms.get(atom[0])][row]+=coef*coef1*parseInt(atom[1]);
        }else if(braL.test(str[i])){
            let nBra=1;
            while(nBra>0&&i+1<str.length){
                if(braR.test(str[i+1])){
                    nBra--;
                    if(nBra>0) atom+=str[i+1];
                }else if(braL.test(str[i+1])){
                    nBra++;
                    atom+=str[i+1];
                }else{
                    atom+=str[i+1];
                }
                i++;
            }
            let coef2="";
            while(/[0-9]/.test(str[i+1])&&i+1<str.length) coef2+=str[1+i++];
            coef2=(coef2==="" ? 1 : parseInt(coef2));
            setAtomCoef(row,atom,atoms,coefMatrix,coef1*coef2);
        }else{
            console.log("undefined:",str[i]);
        }
    }
};

const pivotGaussJordan=(coef)=>{
    let i,j,min=1,tmp;
    let row,maxRow,col=0;
    for(row=0;row<coef.length&&col<coef[0].length-1;row++){
        // 最大pivotを持つ行を探索
        maxRow=coef.reduce((prev,current,i,arr)=>{
            if(i<=row){
                return row;
            }else if(Math.abs(arr[i][col])>Math.abs(arr[prev][col])){
                return i;
            }else{
                return prev;
            }
        },0);
        // 全部の行でcol列がゼロであれば、次の列を探索
        if(coef[maxRow][col]===0){
            col++,row--;
            continue;
        }
        // pivotを最大値に変更
        [coef[row],coef[maxRow]]=[coef[maxRow],coef[row]];
        // pivotの行をpivot値で割る
        tmp=coef[row][col];
        for(i=0;i<coef[0].length;i++) coef[row][i]/=tmp;
        // 他の行に対してpivot行を引き去る
        for(i=0;i<coef.length;i++){
            if(i===row) continue;
            tmp=coef[i][col];
            for(j=0;j<coef[0].length;j++) coef[i][j]-=coef[row][j]*tmp;
        }
    }
    for(i=0;i<coef.length;i++){
        const norm=coef[i].slice(-1)[0];
        for(j=0;j<coef[0].length;j++){
            // 一番右の値で規格化
            if(norm!==0) coef[i][j]/=norm;
            // 解にするために、逆数をとる
            if(coef[i][j]!==0){
                coef[i][j]=1/coef[i][j];
                if((tmp=Math.abs(coef[i][j]))<min&&tmp>1e-10) min=tmp;
            }
        }
    }
    // 小数解の頻度を下げるため、簡易的に全体の最小値を全体にかけて正規化
    if(coef!==1){
        for(i=0;i<coef.length;i++){
            for(j=0;j<coef[0].length;j++){
                coef[i][j]/=min;
            }
        }
    }
};

