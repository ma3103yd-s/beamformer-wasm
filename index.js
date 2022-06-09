import init, {ULA, Signal}  from "./pkg/beamformer_wasm.js";

let source_index = new Set();

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
      break;
    }
    case 2: {
      phi = signal.s_apes(padd);
      phi = phi.map((x) => x);
      phi = phi.map((x) => Math.log(x));
      break;

    }
  }
  return phi; 
}
// Code to run when number of sensors change
function onInput(ula_slider, signal) {
  let val = ula_slider.value;
  document.getElementById("nsens").innerHTML=val;
  console.log(signal);
  signal.clear();
  signal.set_n_sensor(val);
  signal.random_signal();

}
// Function to call when clicking
function onClick(e, PLOT, traces, signal) {
  let names = ["Beamform", "Capon", "Apes"];
  let colors = ["rgb(27, 158, 119)", "rgb(231,41,138)", "rgb(102,166,30)"];
  let index = 0;
  // Find which radio button is checked
  for(let i=0; i<3; i++) {
    if(document.getElementById(names[i]).checked) index = i;
  }
  // Calculate angle of mouse click
  let x =  (e.offsetX-700);
  let y = e.offsetY - 780;
  let theta = Math.atan2(-y, x)*(180/Math.PI)-90;
  let coh = document.getElementById("Coherent").checked;
  let pert = document.getElementById("Pert").value;
  let n_sensor = document.getElementById("nsensor").value;
  signal.set_coh(coh);
  signal.set_pert(pert);
  signal.set_n_sensor(n_sensor);
  let amp = document.getElementById("amp").value;
  signal.add_source(amp, theta);
  source_index.add(index);


  // Sensor padding.
  let padd =  1024;
  let ff = [];
  for (let i = 0; i < padd; i++) {
      ff[i] = -90.0 + (180.0/padd)*i;
  }
  let r_max = 0;
  let temp_trace = [];
  // Loop through the different sources and use corresponding method.
  for(let ind of source_index) {
      let phi = index_to_signal(signal, padd, ind);
      r_max = Math.max(...phi, r_max);
      let trace2 = {
          r: phi,
          theta: ff,
          mode: 'lines',
          name: names[ind],
          marker: {
              color: colors[ind],
              line: {color: "white"},
              opacity: 0.7,
          },
          type: 'scatterpolar',
      };
      temp_trace.push(trace2);

  }

  /* Current Plotly.js version */
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
  // Change so that the signal is always max distance away.
  traces.map((tr) => tr.r = [r_max*1.3]);
  traces.push(trace);

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
  //traces.concat(temp_trace);
  traces = traces.concat(temp_trace);
  Plotly.newPlot(PLOT, traces, layout);
  traces.pop();

}

init()
  .then(() => {
    const PLOT = document.getElementById('plot');
    document.getElementById("Beamform").checked = true;
    document.getElementById("amp").value = 20;
    let ula_slider = document.getElementById("nsensor");
    ula_slider.value = 10;
    document.getElementById("nsens").innerHTML = 10;
    // Define ULA
    let ula = ULA.new(10, 0.5);
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
    ula_slider.addEventListener("change", () => {onInput(ula_slider, signal)});
    PLOT.addEventListener("click", (e) => {onClick(e, PLOT, traces, signal)});  
  });

