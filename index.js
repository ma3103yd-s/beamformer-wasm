import init, {ULA, Signal}  from "./pkg/beamformer_wasm.js";

let last_index = 0;

function coordToAngle(x, y) {
  let theta = Math.atan2(y, x);
  return theta;
}

function addNewSource(e, traces, signal) {

  //let xx = x * (8/344);
  //let yy = y * (8/360);
  //let r = Math.sqrt(xx*xx + yy*yy);
  if (r > 10) r = 10;

  signal.add_source(20.0, theta);
  //traces.push(trace);
  //console.log([x, y]);
}

function index_to_signal(signal, padd, index) {
  let m = signal.n_sensor();
  let phi = [];
   switch (index) {
    case 0: {
      phi = signal.beamform(padd);
      phi = phi.map((x) => x/m/m);
      phi = phi.map((x) => Math.log(x));
      break;
    }
    case 1: {
      phi = signal.capon(padd);
      phi = phi.map((x) => x);
      phi = phi.map((x) => Math.log(x));
    }
  }
  return phi; 
}

function onClick(e, PLOT, traces, signal) {
  let names = ["Beamform", "Capon"];
  let colors = ["rgb(27, 158, 119)", "rgb(231,41,138)"];

  let n_source = signal.n_sources();
  let index = 0;
  for(let i=0; i<2; i++) {
    if(document.getElementById(names[i]).checked) index = i;
  }
  traces.splice(-1, n_source);
  let x =  (e.offsetX-700);
  let y = e.offsetY - 780;
  let theta = Math.atan2(-y, x)*(180/Math.PI)-90;
  let coh = document.getElementById("Coherent").checked;
  let pert = document.getElementById("Pert").value;

  signal.set_coh(coh);
  signal.set_pert(pert);
  let amp = document.getElementById("amp").value;
  signal.add_source(amp, theta);
  console.log(signal);
  //addNewSource(e, traces, signal);

  let padd =  1024;
  let phi1 = [];
  let phi2 = [];
  let ff = [];
  for (let i = 0; i < padd; i++) {
      ff[i] = -90.0 + (180.0/padd)*i;
  }

  phi1 = index_to_signal(signal, padd, index);



  let r_max = Math.max(...phi1,...phi2);
  /* Current Plotly.js version */
  console.log(r_max);
  let trace = {
    r: [1.3*r_max],
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
  traces.map((tr) => tr.r = [r_max*1.3]);
  /*if (traces!=[]) {

  }*/


  
  traces.push(trace);

  let trace1 = {
      r: phi1,
      theta: ff,
      mode: 'lines',
      name: names[index],
      marker: {
          color: colors[index],
          line: {color: "white"},
          opacity: 0.7,
      },
      type: 'scatterpolar',
  };
  if(last_index!=index) {
    phi2 = index_to_signal(signal, padd, last_index);
    let trace2 = {
        r: phi2,
        theta: ff,
        mode: 'lines',
        name: names[last_index],
        marker: {
            color: colors[last_index],
            line: {color: "white"},
            opacity: 0.7,
        },
        type: 'scatterpolar',
    };
    traces.push(trace2);
  }
  let layout = {
    autosize: false,
    width:1400,
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

  last_index = index;
  Plotly.newPlot(PLOT, traces, layout);

}

init()
  .then(() => {
    const PLOT = document.getElementById('plot');
    document.getElementById("Beamform").checked = true;
    document.getElementById("amp").value = 20;

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
      width:1400,
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

