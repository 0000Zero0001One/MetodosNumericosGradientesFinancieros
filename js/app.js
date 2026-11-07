/* ThermoFinance - motor numérico académico
   GC, GCP-Jacobi, Cholesky y autovalores por rotaciones de Jacobi. */
(() => {
"use strict";
const labels = ["Infraestructura","Desarrollo","Seguridad","Capacitación"];
const scenarios = {
 ideal:{name:"Ideal",A:[[10,-1,-1,0],[-1,9,-1,-1],[-1,-1,8,-1],[0,-1,-1,7]],b:[120,100,90,70]},
 stress:{name:"Bajo estrés financiero",A:[[12,-1.2,-1,0],[-1.2,11,-1.1,-1],[-1,-1.1,9.5,-1.2],[0,-1,-1.2,8.5]],b:[135,92,105,64]},
 ill:{name:"Mal condicionado",A:[[1,0.45,1.5,4.5],[0.45,9,4.5,13.5],[1.5,4.5,100,45],[4.5,13.5,45,900]],b:[40,120,300,700]}
};
let activeScenario="ideal", lastResult=null, solutionChart=null, convergenceChart=null;
const $ = id => document.getElementById(id);
const clone = x => JSON.parse(JSON.stringify(x));
const fmt = (x,d=6) => Number.isFinite(x) ? Number(x).toLocaleString("es-BO",{minimumFractionDigits:0,maximumFractionDigits:d}) : "—";
const norm = v => Math.sqrt(v.reduce((s,x)=>s+x*x,0));
const dot = (a,b) => a.reduce((s,x,i)=>s+x*b[i],0);
const matVec = (A,x) => A.map(row=>dot(row,x));
const transpose = A => A[0].map((_,j)=>A.map(row=>row[j]));
const identity = n => Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0));
function renderMatrix(A,b){
  $("matrixBody").innerHTML=A.map((row,i)=>`<tr><th class="matrix-row-label">Ecuación ${i+1}</th>${row.map((v,j)=>`<td><input class="matrix-input a-input" type="number" step="any" data-i="${i}" data-j="${j}" aria-label="A ${i+1}, ${j+1}" value="${v}"></td>`).join("")}<td><input class="matrix-input b-input" type="number" step="any" data-i="${i}" aria-label="b ${i+1}" value="${b[i]}"></td></tr>`).join("");
}
function readSystem(){
  const A=Array.from({length:4},(_,i)=>Array.from({length:4},(_,j)=>Number(document.querySelector(`.a-input[data-i="${i}"][data-j="${j}"]`).value)));
  const b=Array.from({length:4},(_,i)=>Number(document.querySelector(`.b-input[data-i="${i}"]`).value));
  if(A.flat().concat(b).some(v=>!Number.isFinite(v))) throw new Error("Todos los coeficientes deben ser números finitos.");
  return {A,b};
}
function validateSymmetry(A,tol=1e-10){
  for(let i=0;i<A.length;i++)for(let j=i+1;j<A.length;j++)if(Math.abs(A[i][j]-A[j][i])>tol)return false;
  return true;
}
function cholesky(A){
  const n=A.length,L=Array.from({length:n},()=>Array(n).fill(0));
  for(let i=0;i<n;i++)for(let j=0;j<=i;j++){
    let sum=A[i][j];for(let k=0;k<j;k++)sum-=L[i][k]*L[j][k];
    if(i===j){if(!(sum>1e-12))throw new Error("La matriz no es definida positiva: Cholesky encontró un pivote no positivo.");L[i][j]=Math.sqrt(sum);}
    else L[i][j]=sum/L[j][j];
  }
  return L;
}
function solveCholesky(A,b){
  const L=cholesky(A),n=b.length,y=Array(n).fill(0),x=Array(n).fill(0);
  for(let i=0;i<n;i++){let s=b[i];for(let j=0;j<i;j++)s-=L[i][j]*y[j];y[i]=s/L[i][i];}
  for(let i=n-1;i>=0;i--){let s=y[i];for(let j=i+1;j<n;j++)s-=L[j][i]*x[j];x[i]=s/L[i][i];}
  return x;
}
// Autovalores de una matriz simétrica mediante rotaciones de Jacobi.
function symmetricEigenvalues(input){
  const A=clone(input),n=A.length,maxIter=100*n*n;
  for(let iter=0;iter<maxIter;iter++){
    let p=0,q=1,max=0;
    for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(Math.abs(A[i][j])>max){max=Math.abs(A[i][j]);p=i;q=j;}
    if(max<1e-12)break;
    const tau=(A[q][q]-A[p][p])/(2*A[p][q]);
    const t=(tau>=0?1:-1)/(Math.abs(tau)+Math.sqrt(1+tau*tau));
    const c=1/Math.sqrt(1+t*t),s=t*c,app=A[p][p],aqq=A[q][q],apq=A[p][q];
    A[p][p]=app-t*apq;A[q][q]=aqq+t*apq;A[p][q]=0;A[q][p]=0;
    for(let k=0;k<n;k++)if(k!==p&&k!==q){const akp=A[k][p],akq=A[k][q];A[k][p]=A[p][k]=c*akp-s*akq;A[k][q]=A[q][k]=s*akp+c*akq;}
  }
  return A.map((r,i)=>r[i]).sort((a,b)=>a-b);
}
function conditionNumberSPD(A){
  const eig=symmetricEigenvalues(A),min=Math.min(...eig),max=Math.max(...eig);
  if(min<=0)throw new Error("La matriz tiene un autovalor no positivo.");
  return {value:max/min,eigenvalues:eig};
}
function preconditionedCondition(A){
  const d=A.map((r,i)=>r[i]);
  const B=A.map((row,i)=>row.map((v,j)=>v/Math.sqrt(d[i]*d[j])));
  return conditionNumberSPD(B);
}
function residual(A,x,b){const ax=matVec(A,x);return b.map((v,i)=>v-ax[i]);}
function conjugateGradient(A,b,x0,tol,maxIt,preconditioned=false){
  let x=x0.slice(),r=residual(A,x,b),z,p,rz,rr;
  const bnorm=norm(b)||1,history=[norm(r)/bnorm];
  if(preconditioned){z=r.map((v,i)=>v/A[i][i]);p=z.slice();rz=dot(r,z);}
  else{p=r.slice();rr=dot(r,r);}
  let iter=0;
  while(norm(r)/bnorm>tol && iter<maxIt){
    const Ap=matVec(A,p),den=dot(p,Ap);
    if(!(den>0) || !Number.isFinite(den))throw new Error("El algoritmo detectó una dirección no positiva; revise que A sea definida positiva.");
    const alpha=(preconditioned?rz:rr)/den;
    x=x.map((v,i)=>v+alpha*p[i]);r=r.map((v,i)=>v-alpha*Ap[i]);iter++;
    history.push(norm(r)/bnorm);
    if(norm(r)/bnorm<=tol)break;
    if(preconditioned){
      z=r.map((v,i)=>v/A[i][i]);const rzNew=dot(r,z),beta=rzNew/rz;p=z.map((v,i)=>v+beta*p[i]);rz=rzNew;
    }else{
      const rrNew=dot(r,r),beta=rrNew/rr;p=r.map((v,i)=>v+beta*p[i]);rr=rrNew;
    }
  }
  return {x,iterations:iter,residual:norm(r)/bnorm,history,converged:norm(r)/bnorm<=tol};
}
function errorRows(ref,approx){
  return ref.map((v,i)=>{const abs=Math.abs(v-approx[i]);return {abs,rel:v===0?null:abs/Math.abs(v),pct:v===0?null:abs/Math.abs(v)*100};});
}
function setMessage(message,type="info"){
  $("validationMessage").innerHTML=`<div class="alert alert-${type} mb-0">${message}</div>`;
}
function loadScenario(key){
  activeScenario=key;$("scenarioSelect").value=key;renderMatrix(scenarios[key].A,scenarios[key].b);
  $("validationMessage").innerHTML="";lastResult=null;
}
function getX0(){
  if($("x0Mode").value==="zero")return [0,0,0,0];
  const x=$("customX0").value.split(",").map(v=>Number(v.trim()));
  if(x.length!==4||x.some(v=>!Number.isFinite(v)))throw new Error("El vector inicial personalizado debe contener cuatro números separados por coma.");
  return x;
}
function execute(A,b,opts={}){
  if(!validateSymmetry(A))throw new Error("La matriz A no es simétrica. Debe cumplirse A = Aᵀ.");
  cholesky(A); // SPD validation
  const tol=opts.tol??Number($("tolerance").value),maxIt=opts.maxIt??Number($("maxIterations").value);
  if(!(tol>0&&tol<1))throw new Error("La tolerancia debe ser mayor que 0 y menor que 1.");
  if(!(maxIt>=4&&maxIt<=10000))throw new Error("El máximo de iteraciones debe estar entre 4 y 10000.");
  const x0=opts.x0??getX0(),ref=solveCholesky(A,b);
  const gc=conjugateGradient(A,b,x0,tol,maxIt,false),gcp=conjugateGradient(A,b,x0,tol,maxIt,true);
  const cond=conditionNumberSPD(A),pre=preconditionedCondition(A);
  return {A,b,ref,gc,gcp,cond:cond.value,eigen:cond.eigenvalues,preCond:pre.value,errors:errorRows(ref,gcp.x),tol,maxIt};
}
function renderResults(r){
  $("metricCondition").textContent=fmt(r.cond,3);$("metricPreCondition").textContent=fmt(r.preCond,3);
  $("metricIterations").textContent=`${r.gc.iterations} / ${r.gcp.iterations}`;
  $("metricResidual").textContent=r.gcp.residual.toExponential(2);
  $("solutionTable").querySelector("tbody").innerHTML=labels.map((name,i)=>`<tr><th>x${i+1}<small class="d-block text-secondary">${name}</small></th><td>${fmt(r.ref[i],8)}</td><td>${fmt(r.gc.x[i],8)}</td><td>${fmt(r.gcp.x[i],8)}</td></tr>`).join("");
  $("errorTable").querySelector("tbody").innerHTML=r.errors.map((e,i)=>`<tr><th>x${i+1}</th><td>${fmt(e.abs,9)}</td><td>${e.rel===null?"No definido":fmt(e.rel,9)}</td><td>${e.pct===null?"No definido":fmt(e.pct,7)+"%"}</td></tr>`).join("");
  $("diagnosticTable").querySelector("tbody").innerHTML=[
    ["Convergencia",r.gc.converged?"Sí":"No",r.gcp.converged?"Sí":"No"],
    ["Iteraciones",r.gc.iterations,r.gcp.iterations],
    ["Residuo relativo final",r.gc.residual.toExponential(3),r.gcp.residual.toExponential(3)],
    ["Número de condición",fmt(r.cond,6),fmt(r.preCond,6)],
    ["Error absoluto máximo",fmt(Math.max(...errorRows(r.ref,r.gc.x).map(e=>e.abs)),9),fmt(Math.max(...r.errors.map(e=>e.abs)),9)]
  ].map(row=>`<tr>${row.map(v=>`<td>${v}</td>`).join("")}</tr>`).join("");
  renderCharts(r);
}
function renderCharts(r){
  if(!window.Chart)return;
  if(solutionChart)solutionChart.destroy();if(convergenceChart)convergenceChart.destroy();
  solutionChart=new Chart($("solutionChart"),{type:"bar",data:{labels, datasets:[
    {label:"Referencia (Cholesky)",data:r.ref,backgroundColor:"#183f62"},
    {label:"GC",data:r.gc.x,backgroundColor:"#48a4bb"},
    {label:"GCP",data:r.gcp.x,backgroundColor:"#72b99a"}
  ]},options:{responsive:true,plugins:{legend:{position:"bottom"}},scales:{y:{beginAtZero:true,title:{display:true,text:"Miles de Bs"}}}}});
  const maxLen=Math.max(r.gc.history.length,r.gcp.history.length);
  const labelsIt=Array.from({length:maxLen},(_,i)=>i);
  convergenceChart=new Chart($("convergenceChart"),{type:"line",data:{labels:labelsIt,datasets:[
    {label:"GC",data:r.gc.history,borderColor:"#48a4bb",backgroundColor:"#48a4bb",tension:.15,pointRadius:2},
    {label:"GCP - Jacobi",data:r.gcp.history,borderColor:"#183f62",backgroundColor:"#183f62",tension:.15,pointRadius:2}
  ]},options:{responsive:true,plugins:{legend:{position:"bottom"}},scales:{x:{title:{display:true,text:"Iteración"}},y:{type:"logarithmic",title:{display:true,text:"Residuo relativo (escala logarítmica)"},min:1e-16}}}});
}
function runCurrent(){
  try{const {A,b}=readSystem();const r=execute(A,b);lastResult=r;renderResults(r);
    const status=r.gc.converged&&r.gcp.converged?"success":"warning";
    setMessage(`<strong>Validación correcta.</strong> A es simétrica y definida positiva. Cholesky, GC y GCP fueron ejecutados. ${r.gc.converged&&r.gcp.converged?"Ambos métodos alcanzaron la tolerancia.":"Uno o ambos métodos no alcanzaron la tolerancia dentro del máximo de iteraciones."}`,status);
  }catch(e){setMessage(`<strong>No se pudo ejecutar:</strong> ${e.message}`,"danger");}
}
function qaStatus(ok){return ok?'<span class="status-pill status-pass">APROBADO</span>':'<span class="status-pill status-error">FALLIDO</span>';}
function runQA(){
  const results=[];
  // 1. Carga de escenarios
  try{const runs=Object.keys(scenarios).map(k=>execute(scenarios[k].A,scenarios[k].b,{x0:[0,0,0,0],tol:1e-6,maxIt:100}));
    const ok=runs.every(r=>r.gc.converged&&r.gcp.converged);
    results.push({name:"Carga de escenarios",input:"Ideal, estrés y mal condicionado",expected:"Los tres sistemas SPD deben resolverse.",obtained:`${runs.length} escenarios validados; convergencia: ${runs.map(r=>r.gcp.converged?"sí":"no").join(", ")}.`,ok});
  }catch(e){results.push({name:"Carga de escenarios",input:"Ideal, estrés y mal condicionado",expected:"Los tres sistemas SPD deben resolverse.",obtained:e.message,ok:false});}
  // 2. Matriz no simétrica
  const bad=clone(scenarios.ideal.A);bad[0][1]=bad[0][1]+2;
  const ok2=!validateSymmetry(bad);
  results.push({name:"Matriz no válida",input:"Se altera A₁₂ sin modificar A₂₁.",expected:"Detectar falta de simetría y rechazar el sistema.",obtained:ok2?"La matriz no simétrica fue detectada; el sistema se rechazaría antes de resolver.":"No se detectó la matriz inválida.",ok:ok2});
  // 3. Edición dinámica
  try{const before=solveCholesky(scenarios.ideal.A,scenarios.ideal.b);const b2=scenarios.ideal.b.slice();b2[0]+=15;const after=solveCholesky(scenarios.ideal.A,b2);const changed=norm(before.map((v,i)=>v-after[i]))>1e-8;
    results.push({name:"Edición dinámica",input:"Incrementar b₁ en 15 unidades.",expected:"La solución recalculada debe cambiar.",obtained:changed?"La solución cambió tras modificar b y recalcular.":"La solución no cambió.",ok:changed});
  }catch(e){results.push({name:"Edición dinámica",input:"Modificar b₁.",expected:"La solución debe cambiar.",obtained:e.message,ok:false});}
  // 4. PCG vs CG in ill-conditioned
  try{const r=execute(scenarios.ill.A,scenarios.ill.b,{x0:[0,0,0,0],tol:1e-6,maxIt:100});const ok=r.gcp.converged&&r.gc.converged&&r.gcp.iterations<=r.gc.iterations;
    results.push({name:"Comparación de precondicionamiento",input:"Escenario mal condicionado; Jacobi.",expected:"GCP converge en menor o igual número de iteraciones que GC.",obtained:`GC: ${r.gc.iterations} iteraciones; GCP: ${r.gcp.iterations}; κ(A)=${fmt(r.cond,2)}; κ precond.=${fmt(r.preCond,2)}.`,ok});
  }catch(e){results.push({name:"Comparación de precondicionamiento",input:"Escenario mal condicionado.",expected:"GCP debe mejorar o igualar iteraciones.",obtained:e.message,ok:false});}
  $("qaTable").querySelector("tbody").innerHTML=results.map((r,i)=>`<tr><td>${i+1}</td><td>${r.name}</td><td>${r.input}</td><td>${r.expected}</td><td>${r.obtained}</td><td>${qaStatus(r.ok)}</td></tr>`).join("");
  const passed=results.filter(r=>r.ok).length;
  $("qaSummary").innerHTML=`<div class="alert ${passed===4?"alert-success":"alert-warning"} mb-0"><strong>Resultado QA:</strong> ${passed} de 4 pruebas aprobadas.</div>`;
}
$("scenarioSelect").addEventListener("change",e=>loadScenario(e.target.value));
document.querySelectorAll("[data-load-scenario]").forEach(btn=>btn.addEventListener("click",()=>{loadScenario(btn.dataset.loadScenario);$("laboratorio").scrollIntoView({behavior:"smooth"});}));
$("resetScenario").addEventListener("click",()=>loadScenario($("scenarioSelect").value));
$("runButton").addEventListener("click",runCurrent);
$("runQa").addEventListener("click",runQA);
$("x0Mode").addEventListener("change",e=>{$("customX0Wrap").hidden=e.target.value!=="custom";});
$("mentionSelect").addEventListener("change",e=>{
 const notes={
  "IA y Ciencia de Datos":"Énfasis sugerido: preparación de datos, análisis exploratorio y estimación de covarianzas/correlaciones para una cartera de activos. Este sitio conserva el caso presupuestario general; adapta el caso si tu mención exige cartera.",
  "Ingeniería de Sistemas":"Énfasis sugerido: planificación financiera, reglas de negocio, flujo de caja y distribución de recursos entre unidades.",
  "Desarrollo de Sistemas":"Énfasis sugerido: arquitectura web, validaciones matemáticas de entradas y pruebas funcionales del simulador."
 };
 $("mentionNote").textContent=notes[e.target.value]||"Selecciona tu mención para consultar el énfasis adicional recomendado por la consigna.";
});
loadScenario("ideal");
runCurrent();
})();
