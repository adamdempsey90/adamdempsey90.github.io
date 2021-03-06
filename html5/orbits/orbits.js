

function initVars() {
  ecc = Number(eSlider.value);
  semia = Math.pow(10,Number(aSlider.value));
  PlanetMass =  Math.pow(10,Number(mSlider.value));
//  periapseAngle = Number(wSlider.value)
  showNewEccentricity();
  showNewA();
  showNewMass();
  //showNewPeriapse();
  document.getElementById('myButton').value = 'Stop'



  mtot = 1 + PlanetMass;
  ReducedMass = PlanetMass/mtot;
  sfac = PlanetMass/mtot;

  var cos = 1.0;
  var sin = 0.0;
//  var cos = Math.cos(periapseAngle*degToRad);
//  var sin = Math.sin(periapseAngle*degToRad);

  var radTemp = semia*(1-ecc)/(1+PlanetMass); // Distance of peri
  var vradTemp = -Math.sqrt((1+ecc)*NewtG*mtot/(semia*(1-ecc))); // Speed at peri
  vradTemp /= (1+PlanetMass);

 // vx' = vx Cos[x] - vy Sin[x]
 // vy' = vy Cos[x] + vx Sin[x]


  PlanetX =  radTemp*cos;
  PlanetY =  radTemp*sin;
  PlanetVx = -vradTemp*sin;
  PlanetVy = vradTemp*cos;
//  PlanetVy /= (1+PlanetMass);

  SunX = -PlanetX*PlanetMass;
  SunY = -PlanetY*PlanetMass;
  SunVx = -PlanetVx*PlanetMass;
  SunVy = -PlanetVy*PlanetMass;

  chartx = 0;
  chartx2 = 0;
  charty = PlanetX;
  chartyS = SunX;
  charty2 = PlanetVx*VtoKms;
  charty2S = SunVx*VtoKms;
  chartR = SunVx*VtoKms;


  var aubin = setScale();

  drawSun();
  drawPlanet();
  clearTrails();
  drawBar(aubin);
  Moving = true;

  dt = Math.pow(semia*(1-ecc),1.5)/200.;
  dt2 = dt*dt;
  totaltime = 0;
  count = 0;
  timeResolution = Math.pow(semia,1.5)/PointsPerOrbit;

  clearChart();
  updateChart();

  startIntegration();


}
function setScale() {
  var nb;
  var auBins = [.1, .3, 1,3, 10, 30, 100];
  for (var i=0; i<auBins.length; i++) {
    if (semia <= auBins[i] ) {
      nb = auBins[i];
      AuPerPixel = 1.5*2*auBins[i]*(1+ecc)/theCanvas.width; // Sun
      break;
    }
  }
  return nb;
}
function drawBar(nb) {
  trailContext.beginPath();
  trailContext.moveTo(20, 20);
  trailContext.lineTo(nb/AuPerPixel + 20 , 20);
  trailContext.stroke();
  trailContext.font = "15px sans-serif";
  trailContext.fillText(nb + " Au",nb/AuPerPixel + 25,25);

  trailContext.beginPath();
  trailContext.moveTo(20, 40);
  trailContext.lineTo(1./AuPerPixel + 20 , 40);
  trailContext.stroke();
  trailContext.font = "15px sans-serif";
  trailContext.fillText("Earth's Orbit",1/AuPerPixel + 25,45);

}
function buttonChange() {
  if (Moving) {
    Moving=false;
    document.getElementById('myButton').value = 'Start'
  }
  else {
    Moving =true;
    document.getElementById('myButton').value = 'Stop'
    startIntegration();

  }
}


function calcReducedMass() {
  ReducedMass = PlanetMass/(1 + PlanetMass);
}


function drawSun() {

  var pixelX = theCanvas.width/2 + SunX/AuPerPixel;
  var pixelY = theCanvas.height/2 + SunY/AuPerPixel; // Sun

  theContext.clearRect(0, 0, theCanvas.width, theCanvas.height);
  theContext.beginPath();
  theContext.arc(pixelX,pixelY, 20, 0, 2*Math.PI);
  theContext.fillStyle = "yellow"
  theContext.fill();
  trailContext.fillRect(pixelX-0.5, pixelY-0.5, 1, 1);

  pixelX = theCanvas.width/2;
  pixelY = theCanvas.height/2;
  theContext.beginPath();
  theContext.arc(pixelX,pixelY,3,0,2*Math.PI);
  theContext.fillStyle = "black"
  theContext.fill();

  pixelX = theCanvas.width/2 + PlanetX/AuPerPixel;
  pixelY = theCanvas.height/2 + PlanetY/AuPerPixel;

  theContext.beginPath();
  theContext.arc(pixelX,pixelY, PlanetSize, 0, 2*Math.PI);
  theContext.fillStyle = PlanetColor; //"#2E2EFE";
  theContext.fill();
  trailContext.fillRect(pixelX-0.5, pixelY-0.5, 1, 1);

}

function drawPlanet() {

}

function startIntegration() {
  if (Moving) {
    timer = window.setInterval(moveSystemCheck, 1000/fps);
  }
}
function moveSystemCheck() {
  if (Moving) {
    moveSystem();
  }
}
function moveSystem() {



  var sep = [ PlanetX-SunX, PlanetY-SunY];
  var vsep = [PlanetVx - SunVx, PlanetVy - SunVy];
  var accel = computeAccel(sep);


  sep[0] += vsep[0]*dt + .5*accel[0]*dt2;
  sep[1] += vsep[1]*dt + .5*accel[1]*dt2;

  vsep[0] += .5*accel[0]*dt;
  vsep[1] += .5*accel[1]*dt;

  accel = computeAccel(sep);
  vsep[0] += .5*accel[0]*dt;
  vsep[1] += .5*accel[1]*dt;


  PlanetVx = vsep[0]/mtot;
  PlanetVy = vsep[1]/mtot;
  PlanetX = sep[0]/mtot;
  PlanetY = sep[1]/mtot;

  SunVx = -vsep[0]*sfac;
  SunVy = -vsep[1]*sfac;
  SunX = -sep[0]*sfac;
  SunY = -sep[1]*sfac;
  drawSun();
  drawPlanet();

  totaltime += dt;
  count += dt;

  if (count >= timeResolution) {
    chartx = totaltime;
    charty =  PlanetX;//SunVy*VtoKms;
    chartx2 = PlanetX;
    chartyS = SunX;
    charty2 = PlanetVx*VtoKms;
    charty2S = SunVx*VtoKms;
    chartR = SunVx*VtoKms;
    updateChart();
    count = 0;
  }







}

function computeAccel(sep) {
  var rad = Math.sqrt(sep[0]*sep[0] + sep[1]*sep[1]);
  var rad3 = rad*rad*rad;
  var fac = -NewtG*mtot/rad3;
  return [fac*sep[0],fac*sep[1]]
}

function showNewEccentricity() {
  eReadout.innerHTML = ecc;
}
function showNewA() {
  aReadout.innerHTML = semia.toFixed(2) + '&nbsp; Au'
}
function showNewMass() {

   PlanetSize = 5 + 15*(Number(mSlider.value) + 5.5)/(5.5);

  var lbl,num;
  if (Number(mSlider.value) < -4.3) {
    lbl = '&nbsp Earths';
    num = (SuntoEarth*PlanetMass).toFixed(2);
    PlanetColor = 'blue';
  }
  else if (Number(mSlider.value) < -3.5) {
    lbl = '&nbsp Neptunes';
    num = (SuntoNeptune*PlanetMass).toFixed(2);
    PlanetColor = 'green'
  }
  else if (Number(mSlider.value) < -1) {
    lbl = '&nbsp Jupiters';
    num = (SuntoJupiter*PlanetMass).toFixed(2);
    PlanetColor = 'red';
  }
  else {
    lbl = '&nbsp Suns';
    num = PlanetMass.toFixed(2);
    PlanetColor = 'yellow';

  }
  mReadout.innerHTML = num + lbl
}
// function showNewPeriapse() {
//   wReadout.innerHTML = periapseAngle;
// }


function clearTrails() {
  trailContext.clearRect(0, 0, theCanvas.width, theCanvas.height);
}


function clearChart() {
  try {
    chart.destroy();
    chart2.destroy();
    chart3.destroy();

  }
  catch (err) {
  }

  dps = [{x:chartx,y:charty}];
  dps2 = [{x:chartx,y:charty2}];
  dpsS = [{x:chartx,y:chartyS}];
  dps2S = [{x:chartx,y:charty2S}];
  dpsR = [{x:chartx,y:chartR}];
  chart = new CanvasJS.Chart("chartContainer",{
    theme: "theme1",
    backgroundColor: "#FFFFFF",
    title :{
      text: "X Position of Planet (blue) and Sun (red) in Au vs Time (yrs)",
    },
    axisX :{
      title: "Time (yrs)",
    },
    axisY : {
      title: 'X Position (AU)',
    },
    data: [{
      label:"X",
      type: "line",
      dataPoints: dps
    },{
      label:"Y",
      type: "line",
      dataPoints: dpsS,
      color: 'orange'
    }]
  });
  chart2 = new CanvasJS.Chart("chartContainer2",{
    theme: "theme1",
    backgroundColor: "#FFFFFF",
    title :{
      text: "X Velocity of Planet (blue) and Sun (red) in km/s vs Time (yrs)",
    },
    axisX :{
      title: "Time (yrs)",
    },
    axisY : {
      title: 'X Velocity (km/s)',
    },
    data: [{
      label:"V",
      type: "line",
      dataPoints: dps2
    },{
      label:"VS",
      type: "line",
      dataPoints: dps2S,
      color: 'orange'
    }]
  });
  chart3 = new CanvasJS.Chart("chartContainer3",{
    theme: "theme1",
    backgroundColor: "#FFFFFF",
    title :{
      text: "Radial Velocity of Star (km/s) vs Time (yrs)",
    },
    axisX :{
      title: "Time (yrs)",
    },
    axisY : {
      title: 'Radial Velocity (km/s)',
    },
    data: [{
      label:"V",
      type: "line",
      dataPoints: dpsR,
      color: "orange"
    }]
  });

}


function updateChart() {
  // count is number of times loop runs to generate random dataPoints.
    dps.push({
      x: chartx,
      y: charty
    });
    dps2.push({
      x: chartx,
      y: charty2
    });
    dpsS.push({
      x: chartx,
      y: chartyS
    });
    dps2S.push({
      x: chartx,
      y: charty2S
    });
    dpsR.push({
      x: chartx,
      y: chartR

    });


  if (dps.length > dataLength)
  {
    dps.shift();
    dps2.shift();
    dpsS.shift();
    dps2S.shift();
    dpsR.shift();
  }

  chart.render();
  chart2.render();
  chart3.render();

};
