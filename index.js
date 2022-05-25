import init, {ULA, Signal}  from "./pkg/beamformer_wasm.js";

function coordToAngle(x, y) {
  let theta = Math.atan2(y, x);
  return theta;
}

function addNewSource(e, traces, signal) {
  let x =  (e.offsetX-900);
  let y = e.offsetY - 820;
  let theta = Math.atan2(-y, x)*(180/Math.PI)-90;
  let xx = x * (8/344);
  let yy = y * (8/360);
  let r = Math.sqrt(xx*xx + yy*yy);
  if (r > 10) r = 10;

  signal.add_source(20.0, theta);
  let trace = {
    r: [r],
    theta: [theta],
    name: "signal",
    mode: 'marker',
    marker: {
      color: "rgb(217,95,2)", 
      size: 20,
      line: {color: "white"},
      "opacity": 0.7,
    },

    type: "scatterpolar"
  };
  traces.push(trace);
  console.log([x, y]);
}

function onClick(e, PLOT, traces, signal) {
  traces.pop();
  addNewSource(e, traces, signal);
  let padd =  1024;
  let phi = signal.beamform(padd);
  let ff = [];
  for (let i = 0; i < padd; i++) {
      ff[i] = -90.0 + (180.0/padd)*i;
  }
  phi = phi.map((x) => x/400.0);
  phi = phi.map((x) => Math.log(x));
  /* Current Plotly.js version */
  //§console.log(ff);

  let trace1 = {
      r: phi,
      theta: ff,
      mode: 'lines',
      name: 'Beamformer',
      marker: {
          color: "rgb(27, 158, 119)",
          line: {color: "white"},
          opacity: 0.7,
      },
      type: 'scatterpolar',
  };
  let layout = {
    autosize: false,
    width:1800,
    height: 900,
    showlegend: true,
    legend: {
      "orientation": "h",
    },
    polar: {
      sector: [0, 180],
      angularaxis: {
        rotation: 90,
      },
    },
    automargin: true,
  };
  traces.push(trace1);
  Plotly.newPlot(PLOT, traces, layout);

}

init()
  .then(() => {
    const PLOT = document.getElementById('plot');
    const ula = ULA.new(20, 0.5);
    let traces = [];
    let signal = Signal.new(ula, 64, [], []);
    let first_trace = {
      r: [],
      theta: [],
      type:  "scatterpolar",
    }
    let layout = {
      autosize: false,
      width:1800,
      height: 900,
      showlegend: true,
      legend: {
        "orientation": "h",
      },
      polar: {
        sector: [0, 180],
        angularaxis: {
          rotation: 90,
        }
      },
      automargin: true,
    };
    Plotly.newPlot(PLOT, [first_trace], layout);
    PLOT.addEventListener("click", (e) => {onClick(e, PLOT, traces, signal)});
    //let doa = [25.0, 36.0];
    /*let padd = 1024
    //let sigAmp = [10.0, 20.0];

    //let test_signal = Signal.new(ula, 64, sigAmp, doa).with_random_signal();
    let phi = signal.beamform(1024);

    let ff = [];
    for (let i = 0; i < padd; i++) {
        ff[i] = -90.0 + (180.0/padd)*i;
    }
    phi = phi.map((x) => x/400.0);
    phi = phi.map((x) => Math.log(x));
  
    let trace1 = {
        r: phi,
        theta: ff,
        mode: 'lines',
        name: 'Beamformer',
        marker: {
            color: "rgb(27, 158, 119)",
            line: {color: "white"},
            opacity: 0.7,
        },
        type: 'scatterpolargl',
    };
    let layout = {
      autosize: false,
      width:1000,
      height: 600,
      showlegend: true,
      legend: {
        "orientation": "h",
      },
      polar: {
        angularaxis: {
          rotation: 90,
        }
      },
      automargin: true,
    };
    traces.push(trace1);
    Plotly.newPlot(PLOT, traces, layout);
    */
   
  });

