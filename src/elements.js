/* BlueMolecule — element data.
 * Covalent radii: Cordero et al., Dalton Trans. 2008 (low-spin values for Mn/Fe/Co).
 * Van der Waals radii: Bondi 1964 with common extensions; 2.0 Å where unknown.
 * Colours: the widely used CPK/Jmol-style palette.
 * Masses: IUPAC standard atomic weights (abridged).
 */
(function (BM) {
  'use strict';
  const SYMBOLS = ['X',
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'];

  const NAMES = ['Dummy',
    'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron', 'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon',
    'Sodium', 'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur', 'Chlorine', 'Argon',
    'Potassium', 'Calcium', 'Scandium', 'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron', 'Cobalt', 'Nickel',
    'Copper', 'Zinc', 'Gallium', 'Germanium', 'Arsenic', 'Selenium', 'Bromine', 'Krypton',
    'Rubidium', 'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum', 'Technetium', 'Ruthenium', 'Rhodium',
    'Palladium', 'Silver', 'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium', 'Iodine', 'Xenon',
    'Caesium', 'Barium', 'Lanthanum', 'Cerium', 'Praseodymium', 'Neodymium', 'Promethium', 'Samarium', 'Europium',
    'Gadolinium', 'Terbium', 'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium', 'Lutetium',
    'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium', 'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury', 'Thallium',
    'Lead', 'Bismuth', 'Polonium', 'Astatine', 'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium', 'Protactinium',
    'Uranium', 'Neptunium', 'Plutonium', 'Americium', 'Curium', 'Berkelium', 'Californium', 'Einsteinium', 'Fermium',
    'Mendelevium', 'Nobelium', 'Lawrencium', 'Rutherfordium', 'Dubnium', 'Seaborgium', 'Bohrium', 'Hassium',
    'Meitnerium', 'Darmstadtium', 'Roentgenium', 'Copernicium', 'Nihonium', 'Flerovium', 'Moscovium', 'Livermorium',
    'Tennessine', 'Oganesson'];

  const COV = [0.30,
    0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
    1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
    2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16,
    2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40,
    2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87,
    1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50,
    2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69];

  const VDW_KNOWN = {
    1: 1.20, 2: 1.40, 3: 1.82, 4: 1.53, 5: 1.92, 6: 1.70, 7: 1.55, 8: 1.52, 9: 1.47, 10: 1.54,
    11: 2.27, 12: 1.73, 13: 1.84, 14: 2.10, 15: 1.80, 16: 1.80, 17: 1.75, 18: 1.88, 19: 2.75, 20: 2.31,
    28: 1.63, 29: 1.40, 30: 1.39, 31: 1.87, 32: 2.11, 33: 1.85, 34: 1.90, 35: 1.85, 36: 2.02,
    37: 3.03, 38: 2.49, 46: 1.63, 47: 1.72, 48: 1.58, 49: 1.93, 50: 2.17, 51: 2.06, 52: 2.06, 53: 1.98, 54: 2.16,
    55: 3.43, 56: 2.68, 78: 1.75, 79: 1.66, 80: 1.55, 81: 1.96, 82: 2.02, 83: 2.07, 84: 1.97, 85: 2.02, 86: 2.20,
    87: 3.48, 88: 2.83, 92: 1.86
  };

  const COLORS_HEX = {
    1: 'FFFFFF', 2: 'D9FFFF', 3: 'CC80FF', 4: 'C2FF00', 5: 'FFB5B5', 6: '909090', 7: '3050F8', 8: 'FF0D0D', 9: '90E050', 10: 'B3E3F5',
    11: 'AB5CF2', 12: '8AFF00', 13: 'BFA6A6', 14: 'F0C8A0', 15: 'FF8000', 16: 'FFFF30', 17: '1FF01F', 18: '80D1E3',
    19: '8F40D4', 20: '3DFF00', 21: 'E6E6E6', 22: 'BFC2C7', 23: 'A6A6AB', 24: '8A99C7', 25: '9C7AC7', 26: 'E06633', 27: 'F090A0',
    28: '50D050', 29: 'C88033', 30: '7D80B0', 31: 'C28F8F', 32: '668F8F', 33: 'BD80E3', 34: 'FFA100', 35: 'A62929', 36: '5CB8D1',
    37: '702EB0', 38: '00FF00', 39: '94FFFF', 40: '94E0E0', 41: '73C2C9', 42: '54B5B5', 43: '3B9E9E', 44: '248F8F', 45: '0A7D8C',
    46: '006985', 47: 'C0C0C0', 48: 'FFD98F', 49: 'A67573', 50: '668080', 51: '9E63B5', 52: 'D47A00', 53: '940094', 54: '429EB0',
    55: '57178F', 56: '00C900', 57: '70D4FF', 58: 'FFFFC7', 59: 'D9FFC7', 60: 'C7FFC7', 61: 'A3FFC7', 62: '8FFFC7', 63: '61FFC7',
    64: '45FFC7', 65: '30FFC7', 66: '1FFFC7', 67: '00FF9C', 68: '00E675', 69: '00D452', 70: '00BF38', 71: '00AB24',
    72: '4DC2FF', 73: '4DA6FF', 74: '2194D6', 75: '267DAB', 76: '266696', 77: '175487', 78: 'D0D0E0', 79: 'FFD123', 80: 'B8B8D0',
    81: 'A6544D', 82: '575961', 83: '9E4FB5', 84: 'AB5C00', 85: '754F45', 86: '428296', 87: '420066', 88: '007D00', 89: '70ABFA',
    90: '00BAFF', 91: '00A1FF', 92: '008FFF', 93: '0080FF', 94: '006BFF', 95: '545CF2', 96: '785CE3', 97: '8A4FE3', 98: 'A136D4',
    99: 'B31FD4', 100: 'B31FBA', 101: 'B30DA6', 102: 'BD0D87', 103: 'C70066', 104: 'CC0059', 105: 'D1004F', 106: 'D90045',
    107: 'E00038', 108: 'E6002E', 109: 'EB0026'
  };

  const MASS = [0,
    1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180,
    22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948,
    39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38,
    69.723, 72.630, 74.922, 78.971, 79.904, 83.798,
    85.468, 87.62, 88.906, 91.224, 92.906, 95.95, 98, 101.07, 102.91, 106.42, 107.87, 112.41,
    114.82, 118.71, 121.76, 127.60, 126.90, 131.29,
    132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 145, 150.36, 151.96, 157.25, 158.93, 162.50, 164.93,
    167.26, 168.93, 173.05, 174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59,
    204.38, 207.2, 208.98, 209, 210, 222,
    223, 226, 227, 232.04, 231.04, 238.03, 237, 244, 243, 247, 247, 251, 252, 257, 258, 259, 266,
    267, 268, 269, 270, 277, 278, 281, 282, 285, 286, 289, 290, 293, 294, 294];

  const N = SYMBOLS.length;
  const covalent = new Float32Array(N);
  const vdw = new Float32Array(N);
  const color = [];
  const mass = new Float64Array(N);
  const lookup = Object.create(null);
  for (let z = 0; z < N; z++) {
    covalent[z] = z < COV.length ? COV[z] : 1.6;
    vdw[z] = VDW_KNOWN[z] || (z === 0 ? 1.0 : 2.0);
    const hex = COLORS_HEX[z] || (z === 0 ? 'FF66CC' : 'FF1493');
    color.push([parseInt(hex.slice(0, 2), 16) / 255, parseInt(hex.slice(2, 4), 16) / 255, parseInt(hex.slice(4, 6), 16) / 255]);
    mass[z] = z < MASS.length ? MASS[z] : 0;
    lookup[SYMBOLS[z].toLowerCase()] = z;
  }
  // Dummy / ghost labels used by quantum chemistry programs.
  ['x', 'xx', 'bq', 'gh', 'du', 'da', 'q'].forEach(s => { lookup[s] = 0; });

  /** Atomic number from a symbol or label ("C", "CL", "Cl12", "C-CA", "6", "Fe(Iso=56)"). Returns -1 if unknown. */
  function toZ(label) {
    if (label === null || label === undefined) return -1;
    let s = String(label).trim();
    if (!s) return -1;
    if (/^\d+$/.test(s)) { const z = parseInt(s, 10); return z < N ? z : -1; }
    const m = s.match(/^([A-Za-z])([A-Za-z]?)/);
    if (!m) return -1;
    const one = m[1].toLowerCase();
    if (m[2]) {
      const two = (m[1] + m[2]).toLowerCase();
      const rest = s.slice(2);
      const secondIsLower = m[2] === m[2].toLowerCase();
      // Mixed case ("Cl", "Fe") or an all-caps two-letter symbol with nothing else alphabetic after it ("CL", "CL1").
      if (two in lookup && (secondIsLower || !/^[A-Za-z]/.test(rest))) {
        // Avoid reading "HC1"/"CA12"-style force-field labels as two-letter elements only when lower-case; caps handled above.
        return lookup[two];
      }
    }
    if (one in lookup) return lookup[one];
    return -1;
  }

  function symbol(z) { return z >= 0 && z < N ? SYMBOLS[z] : '?'; }
  function name(z) { return z >= 0 && z < N ? NAMES[z] : 'Unknown'; }

  BM.el = { SYMBOLS, NAMES, covalent, vdw, color, mass, toZ, symbol, name, count: N };
})(typeof window !== 'undefined' ? (window.BM = window.BM || {}) : (globalThis.BM = globalThis.BM || {}));
