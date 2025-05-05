/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Atomic database for element properties
//!
//! This module provides atomic data such as element symbols, atomic weights,
//! and other physical properties needed for FEFF calculations.

/// Provides element symbols for atomic numbers
pub fn element_symbol(atomic_number: i32) -> Option<&'static str> {
    match atomic_number {
        1 => Some("H"),
        2 => Some("He"),
        3 => Some("Li"),
        4 => Some("Be"),
        5 => Some("B"),
        6 => Some("C"),
        7 => Some("N"),
        8 => Some("O"),
        9 => Some("F"),
        10 => Some("Ne"),
        11 => Some("Na"),
        12 => Some("Mg"),
        13 => Some("Al"),
        14 => Some("Si"),
        15 => Some("P"),
        16 => Some("S"),
        17 => Some("Cl"),
        18 => Some("Ar"),
        19 => Some("K"),
        20 => Some("Ca"),
        21 => Some("Sc"),
        22 => Some("Ti"),
        23 => Some("V"),
        24 => Some("Cr"),
        25 => Some("Mn"),
        26 => Some("Fe"),
        27 => Some("Co"),
        28 => Some("Ni"),
        29 => Some("Cu"),
        30 => Some("Zn"),
        31 => Some("Ga"),
        32 => Some("Ge"),
        33 => Some("As"),
        34 => Some("Se"),
        35 => Some("Br"),
        36 => Some("Kr"),
        37 => Some("Rb"),
        38 => Some("Sr"),
        39 => Some("Y"),
        40 => Some("Zr"),
        41 => Some("Nb"),
        42 => Some("Mo"),
        43 => Some("Tc"),
        44 => Some("Ru"),
        45 => Some("Rh"),
        46 => Some("Pd"),
        47 => Some("Ag"),
        48 => Some("Cd"),
        49 => Some("In"),
        50 => Some("Sn"),
        51 => Some("Sb"),
        52 => Some("Te"),
        53 => Some("I"),
        54 => Some("Xe"),
        55 => Some("Cs"),
        56 => Some("Ba"),
        57 => Some("La"),
        58 => Some("Ce"),
        59 => Some("Pr"),
        60 => Some("Nd"),
        61 => Some("Pm"),
        62 => Some("Sm"),
        63 => Some("Eu"),
        64 => Some("Gd"),
        65 => Some("Tb"),
        66 => Some("Dy"),
        67 => Some("Ho"),
        68 => Some("Er"),
        69 => Some("Tm"),
        70 => Some("Yb"),
        71 => Some("Lu"),
        72 => Some("Hf"),
        73 => Some("Ta"),
        74 => Some("W"),
        75 => Some("Re"),
        76 => Some("Os"),
        77 => Some("Ir"),
        78 => Some("Pt"),
        79 => Some("Au"),
        80 => Some("Hg"),
        81 => Some("Tl"),
        82 => Some("Pb"),
        83 => Some("Bi"),
        84 => Some("Po"),
        85 => Some("At"),
        86 => Some("Rn"),
        87 => Some("Fr"),
        88 => Some("Ra"),
        89 => Some("Ac"),
        90 => Some("Th"),
        91 => Some("Pa"),
        92 => Some("U"),
        93 => Some("Np"),
        94 => Some("Pu"),
        95 => Some("Am"),
        96 => Some("Cm"),
        97 => Some("Bk"),
        98 => Some("Cf"),
        99 => Some("Es"),
        100 => Some("Fm"),
        101 => Some("Md"),
        102 => Some("No"),
        103 => Some("Lr"),
        104 => Some("Rf"),
        105 => Some("Db"),
        106 => Some("Sg"),
        107 => Some("Bh"),
        108 => Some("Hs"),
        109 => Some("Mt"),
        110 => Some("Ds"),
        111 => Some("Rg"),
        112 => Some("Cn"),
        113 => Some("Nh"),
        114 => Some("Fl"),
        115 => Some("Mc"),
        116 => Some("Lv"),
        117 => Some("Ts"),
        118 => Some("Og"),
        _ => None,
    }
}

/// Returns the atomic weight in atomic mass units (amu)
///
/// Values are based on the relative atomic masses from IUPAC 2013
pub fn atomic_weight(atomic_number: i32) -> Option<f64> {
    match atomic_number {
        1 => Some(1.008),
        2 => Some(4.0026),
        3 => Some(6.94),
        4 => Some(9.0122),
        5 => Some(10.81),
        6 => Some(12.011),
        7 => Some(14.007),
        8 => Some(15.999),
        9 => Some(18.998),
        10 => Some(20.180),
        11 => Some(22.990),
        12 => Some(24.305),
        13 => Some(26.982),
        14 => Some(28.085),
        15 => Some(30.974),
        16 => Some(32.06),
        17 => Some(35.45),
        18 => Some(39.95),
        19 => Some(39.098),
        20 => Some(40.078),
        21 => Some(44.956),
        22 => Some(47.867),
        23 => Some(50.942),
        24 => Some(51.996),
        25 => Some(54.938),
        26 => Some(55.845),
        27 => Some(58.933),
        28 => Some(58.693),
        29 => Some(63.546),
        30 => Some(65.38),
        31 => Some(69.723),
        32 => Some(72.630),
        33 => Some(74.922),
        34 => Some(78.971),
        35 => Some(79.904),
        36 => Some(83.798),
        37 => Some(85.468),
        38 => Some(87.62),
        39 => Some(88.906),
        40 => Some(91.224),
        41 => Some(92.906),
        42 => Some(95.95),
        43 => Some(98.0), // Technetium - no stable isotopes, using approx.
        44 => Some(101.07),
        45 => Some(102.91),
        46 => Some(106.42),
        47 => Some(107.87),
        48 => Some(112.41),
        49 => Some(114.82),
        50 => Some(118.71),
        51 => Some(121.76),
        52 => Some(127.60),
        53 => Some(126.90),
        54 => Some(131.29),
        55 => Some(132.91),
        56 => Some(137.33),
        57 => Some(138.91),
        58 => Some(140.12),
        59 => Some(140.91),
        60 => Some(144.24),
        61 => Some(145.0), // Promethium - no stable isotopes, using approx.
        62 => Some(150.36),
        63 => Some(151.96),
        64 => Some(157.25),
        65 => Some(158.93),
        66 => Some(162.50),
        67 => Some(164.93),
        68 => Some(167.26),
        69 => Some(168.93),
        70 => Some(173.05),
        71 => Some(174.97),
        72 => Some(178.49),
        73 => Some(180.95),
        74 => Some(183.84),
        75 => Some(186.21),
        76 => Some(190.23),
        77 => Some(192.22),
        78 => Some(195.08),
        79 => Some(196.97),
        80 => Some(200.59),
        81 => Some(204.38),
        82 => Some(207.2),
        83 => Some(208.98),
        84 => Some(209.0), // Polonium - approximated
        85 => Some(210.0), // Astatine - approximated
        86 => Some(222.0), // Radon - approximated
        87 => Some(223.0), // Francium - approximated
        88 => Some(226.0), // Radium - approximated
        89 => Some(227.0), // Actinium - approximated
        90 => Some(232.04),
        91 => Some(231.04),
        92 => Some(238.03),
        93 => Some(237.0),  // Neptunium - approximated
        94 => Some(244.0),  // Plutonium - approximated
        95 => Some(243.0),  // Americium - approximated
        96 => Some(247.0),  // Curium - approximated
        97 => Some(247.0),  // Berkelium - approximated
        98 => Some(251.0),  // Californium - approximated
        99 => Some(252.0),  // Einsteinium - approximated
        100 => Some(257.0), // Fermium - approximated
        101 => Some(258.0), // Mendelevium - approximated
        102 => Some(259.0), // Nobelium - approximated
        103 => Some(266.0), // Lawrencium - approximated
        104 => Some(267.0), // Rutherfordium - approximated
        105 => Some(268.0), // Dubnium - approximated
        106 => Some(269.0), // Seaborgium - approximated
        107 => Some(270.0), // Bohrium - approximated
        108 => Some(277.0), // Hassium - approximated
        109 => Some(278.0), // Meitnerium - approximated
        110 => Some(281.0), // Darmstadtium - approximated
        111 => Some(282.0), // Roentgenium - approximated
        112 => Some(285.0), // Copernicium - approximated
        113 => Some(286.0), // Nihonium - approximated
        114 => Some(289.0), // Flerovium - approximated
        115 => Some(290.0), // Moscovium - approximated
        116 => Some(293.0), // Livermorium - approximated
        117 => Some(294.0), // Tennessine - approximated
        118 => Some(294.0), // Oganesson - approximated
        _ => None,
    }
}

/// Returns the covalent radius in Angstroms
pub fn covalent_radius(atomic_number: i32) -> Option<f64> {
    match atomic_number {
        1 => Some(0.31),  // H
        2 => Some(0.28),  // He
        3 => Some(1.28),  // Li
        4 => Some(0.96),  // Be
        5 => Some(0.84),  // B
        6 => Some(0.76),  // C
        7 => Some(0.71),  // N
        8 => Some(0.66),  // O
        9 => Some(0.57),  // F
        10 => Some(0.58), // Ne
        11 => Some(1.66), // Na
        12 => Some(1.41), // Mg
        13 => Some(1.21), // Al
        14 => Some(1.11), // Si
        15 => Some(1.07), // P
        16 => Some(1.05), // S
        17 => Some(1.02), // Cl
        18 => Some(1.06), // Ar
        19 => Some(2.03), // K
        20 => Some(1.76), // Ca
        21 => Some(1.70), // Sc
        22 => Some(1.60), // Ti
        23 => Some(1.53), // V
        24 => Some(1.39), // Cr
        25 => Some(1.39), // Mn
        26 => Some(1.32), // Fe
        27 => Some(1.26), // Co
        28 => Some(1.24), // Ni
        29 => Some(1.32), // Cu
        30 => Some(1.22), // Zn
        31 => Some(1.22), // Ga
        32 => Some(1.20), // Ge
        33 => Some(1.19), // As
        34 => Some(1.20), // Se
        35 => Some(1.20), // Br
        36 => Some(1.16), // Kr
        37 => Some(2.20), // Rb
        38 => Some(1.95), // Sr
        39 => Some(1.90), // Y
        40 => Some(1.75), // Zr
        41 => Some(1.64), // Nb
        42 => Some(1.54), // Mo
        43 => Some(1.47), // Tc
        44 => Some(1.46), // Ru
        45 => Some(1.42), // Rh
        46 => Some(1.39), // Pd
        47 => Some(1.45), // Ag
        48 => Some(1.44), // Cd
        49 => Some(1.42), // In
        50 => Some(1.39), // Sn
        51 => Some(1.39), // Sb
        52 => Some(1.38), // Te
        53 => Some(1.39), // I
        54 => Some(1.40), // Xe
        55 => Some(2.44), // Cs
        56 => Some(2.15), // Ba
        57 => Some(2.07), // La
        58 => Some(2.04), // Ce
        59 => Some(2.03), // Pr
        60 => Some(2.01), // Nd
        61 => Some(1.99), // Pm
        62 => Some(1.98), // Sm
        63 => Some(1.98), // Eu
        64 => Some(1.96), // Gd
        65 => Some(1.94), // Tb
        66 => Some(1.92), // Dy
        67 => Some(1.92), // Ho
        68 => Some(1.89), // Er
        69 => Some(1.90), // Tm
        70 => Some(1.87), // Yb
        71 => Some(1.87), // Lu
        72 => Some(1.75), // Hf
        73 => Some(1.70), // Ta
        74 => Some(1.62), // W
        75 => Some(1.51), // Re
        76 => Some(1.44), // Os
        77 => Some(1.41), // Ir
        78 => Some(1.36), // Pt
        79 => Some(1.36), // Au
        80 => Some(1.32), // Hg
        81 => Some(1.45), // Tl
        82 => Some(1.46), // Pb
        83 => Some(1.48), // Bi
        84 => Some(1.40), // Po
        85 => Some(1.50), // At
        86 => Some(1.50), // Rn
        87 => Some(2.60), // Fr
        88 => Some(2.21), // Ra
        89 => Some(2.15), // Ac
        90 => Some(2.06), // Th
        91 => Some(2.00), // Pa
        92 => Some(1.96), // U
        93 => Some(1.90), // Np
        94 => Some(1.87), // Pu
        95 => Some(1.80), // Am
        96 => Some(1.69), // Cm
        // Approximated values for the rest
        97..=118 => Some(1.65),
        _ => None,
    }
}

/// Returns the K-edge energy in eV
pub fn k_edge_energy(atomic_number: i32) -> Option<f64> {
    match atomic_number {
        1 => None, // H - K-edge not typically measured
        2 => None, // He - K-edge not typically measured
        3 => Some(54.7),
        4 => Some(111.5),
        5 => Some(188.0),
        6 => Some(284.2),
        7 => Some(409.9),
        8 => Some(543.1),
        9 => Some(696.7),
        10 => Some(870.2),
        11 => Some(1070.8),
        12 => Some(1303.0),
        13 => Some(1559.6),
        14 => Some(1839.0),
        15 => Some(2145.5),
        16 => Some(2472.0),
        17 => Some(2822.4),
        18 => Some(3205.9),
        19 => Some(3608.4),
        20 => Some(4038.5),
        21 => Some(4492.0),
        22 => Some(4966.0),
        23 => Some(5465.0),
        24 => Some(5989.0),
        25 => Some(6539.0),
        26 => Some(7112.0),
        27 => Some(7709.0),
        28 => Some(8333.0),
        29 => Some(8979.0),
        30 => Some(9659.0),
        31 => Some(10367.0),
        32 => Some(11103.0),
        33 => Some(11867.0),
        34 => Some(12658.0),
        35 => Some(13474.0),
        36 => Some(14326.0),
        37 => Some(15200.0),
        38 => Some(16105.0),
        39 => Some(17038.0),
        40 => Some(17998.0),
        41 => Some(18986.0),
        42 => Some(20000.0),
        43 => Some(21044.0),
        44 => Some(22117.0),
        45 => Some(23220.0),
        46 => Some(24350.0),
        47 => Some(25514.0),
        48 => Some(26711.0),
        49 => Some(27940.0),
        50 => Some(29200.0),
        51 => Some(30491.0),
        52 => Some(31814.0),
        53 => Some(33169.0),
        54 => Some(34561.0),
        55 => Some(35985.0),
        56 => Some(37441.0),
        57 => Some(38925.0),
        58 => Some(40443.0),
        59 => Some(41991.0),
        60 => Some(43569.0),
        61 => Some(45184.0),
        62 => Some(46834.0),
        63 => Some(48519.0),
        64 => Some(50239.0),
        65 => Some(51996.0),
        66 => Some(53789.0),
        67 => Some(55618.0),
        68 => Some(57486.0),
        69 => Some(59390.0),
        70 => Some(61332.0),
        71 => Some(63314.0),
        72 => Some(65351.0),
        73 => Some(67416.0),
        74 => Some(69525.0),
        75 => Some(71676.0),
        76 => Some(73871.0),
        77 => Some(76111.0),
        78 => Some(78395.0),
        79 => Some(80725.0),
        80 => Some(83102.0),
        81 => Some(85530.0),
        82 => Some(88005.0),
        83 => Some(90526.0),
        84 => Some(93105.0),
        85 => Some(95730.0),
        86 => Some(98404.0),
        87 => Some(101137.0),
        88 => Some(103922.0),
        89 => Some(106755.0),
        90 => Some(109651.0),
        91 => Some(112601.0),
        92 => Some(115606.0),
        93 => Some(118669.0),
        94 => Some(121791.0),
        95 => Some(124982.0),
        96 => Some(128241.0),
        97 => Some(131575.0),
        98 => Some(134987.0),
        // Approximate for remaining elements
        99..=118 => atomic_number.checked_mul(1400).map(|v| v as f64),
        _ => None,
    }
}

/// Returns the L1-edge energy in eV
pub fn l1_edge_energy(atomic_number: i32) -> Option<f64> {
    // Only provide data for elements that have well-defined L1 edges (Z >= 21)
    match atomic_number {
        21 => Some(498.0),  // Sc
        22 => Some(563.0),  // Ti
        23 => Some(628.0),  // V
        24 => Some(695.0),  // Cr
        25 => Some(769.0),  // Mn
        26 => Some(846.0),  // Fe
        27 => Some(926.0),  // Co
        28 => Some(1008.0), // Ni
        29 => Some(1096.0), // Cu
        30 => Some(1193.0), // Zn
        31 => Some(1297.0), // Ga
        32 => Some(1414.0), // Ge
        33 => Some(1527.0), // As
        34 => Some(1652.0), // Se
        35 => Some(1782.0), // Br
        36 => Some(1921.0), // Kr
        37 => Some(2065.0), // Rb
        38 => Some(2216.0), // Sr
        39 => Some(2373.0), // Y
        40 => Some(2532.0), // Zr
        41 => Some(2698.0), // Nb
        42 => Some(2866.0), // Mo
        43 => Some(3043.0), // Tc
        44 => Some(3224.0), // Ru
        45 => Some(3412.0), // Rh
        46 => Some(3604.0), // Pd
        47 => Some(3806.0), // Ag
        48 => Some(4018.0), // Cd
        49 => Some(4238.0), // In
        50 => Some(4465.0), // Sn
        51 => Some(4698.0), // Sb
        52 => Some(4939.0), // Te
        53 => Some(5188.0), // I
        54 => Some(5453.0), // Xe
        55 => Some(5714.0), // Cs
        56 => Some(5989.0), // Ba
        57 => Some(6266.0), // La
        58 => Some(6548.0), // Ce
        // Remaining elements (approximated)
        59..=118 => Some((atomic_number as f64 - 36.0) * 290.0),
        _ => None,
    }
}

/// Returns the atomic number for an element symbol
///
/// This function is case-insensitive and will handle both "Fe" and "FE"
pub fn atomic_number_from_symbol(symbol: &str) -> Option<i32> {
    // Convert to title case for consistent lookup
    let symbol = symbol.to_lowercase();
    let symbol = match symbol.len() {
        1 => symbol.to_uppercase(),
        2..=3 => {
            let mut chars = symbol.chars();
            let first = chars.next().unwrap().to_uppercase().to_string();
            let rest: String = chars.collect();
            format!("{}{}", first, rest)
        }
        _ => return None,
    };

    // Handle special case for two-character symbols with trailing whitespace
    let symbol = symbol.trim();

    match symbol {
        "H" => Some(1),
        "He" => Some(2),
        "Li" => Some(3),
        "Be" => Some(4),
        "B" => Some(5),
        "C" => Some(6),
        "N" => Some(7),
        "O" => Some(8),
        "F" => Some(9),
        "Ne" => Some(10),
        "Na" => Some(11),
        "Mg" => Some(12),
        "Al" => Some(13),
        "Si" => Some(14),
        "P" => Some(15),
        "S" => Some(16),
        "Cl" => Some(17),
        "Ar" => Some(18),
        "K" => Some(19),
        "Ca" => Some(20),
        "Sc" => Some(21),
        "Ti" => Some(22),
        "V" => Some(23),
        "Cr" => Some(24),
        "Mn" => Some(25),
        "Fe" => Some(26),
        "Co" => Some(27),
        "Ni" => Some(28),
        "Cu" => Some(29),
        "Zn" => Some(30),
        "Ga" => Some(31),
        "Ge" => Some(32),
        "As" => Some(33),
        "Se" => Some(34),
        "Br" => Some(35),
        "Kr" => Some(36),
        "Rb" => Some(37),
        "Sr" => Some(38),
        "Y" => Some(39),
        "Zr" => Some(40),
        "Nb" => Some(41),
        "Mo" => Some(42),
        "Tc" => Some(43),
        "Ru" => Some(44),
        "Rh" => Some(45),
        "Pd" => Some(46),
        "Ag" => Some(47),
        "Cd" => Some(48),
        "In" => Some(49),
        "Sn" => Some(50),
        "Sb" => Some(51),
        "Te" => Some(52),
        "I" => Some(53),
        "Xe" => Some(54),
        "Cs" => Some(55),
        "Ba" => Some(56),
        "La" => Some(57),
        "Ce" => Some(58),
        "Pr" => Some(59),
        "Nd" => Some(60),
        "Pm" => Some(61),
        "Sm" => Some(62),
        "Eu" => Some(63),
        "Gd" => Some(64),
        "Tb" => Some(65),
        "Dy" => Some(66),
        "Ho" => Some(67),
        "Er" => Some(68),
        "Tm" => Some(69),
        "Yb" => Some(70),
        "Lu" => Some(71),
        "Hf" => Some(72),
        "Ta" => Some(73),
        "W" => Some(74),
        "Re" => Some(75),
        "Os" => Some(76),
        "Ir" => Some(77),
        "Pt" => Some(78),
        "Au" => Some(79),
        "Hg" => Some(80),
        "Tl" => Some(81),
        "Pb" => Some(82),
        "Bi" => Some(83),
        "Po" => Some(84),
        "At" => Some(85),
        "Rn" => Some(86),
        "Fr" => Some(87),
        "Ra" => Some(88),
        "Ac" => Some(89),
        "Th" => Some(90),
        "Pa" => Some(91),
        "U" => Some(92),
        "Np" => Some(93),
        "Pu" => Some(94),
        "Am" => Some(95),
        "Cm" => Some(96),
        "Bk" => Some(97),
        "Cf" => Some(98),
        "Es" => Some(99),
        "Fm" => Some(100),
        "Md" => Some(101),
        "No" => Some(102),
        "Lr" => Some(103),
        "Rf" => Some(104),
        "Db" => Some(105),
        "Sg" => Some(106),
        "Bh" => Some(107),
        "Hs" => Some(108),
        "Mt" => Some(109),
        "Ds" => Some(110),
        "Rg" => Some(111),
        "Cn" => Some(112),
        "Nh" => Some(113),
        "Fl" => Some(114),
        "Mc" => Some(115),
        "Lv" => Some(116),
        "Ts" => Some(117),
        "Og" => Some(118),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_element_symbol() {
        assert_eq!(element_symbol(1), Some("H"));
        assert_eq!(element_symbol(6), Some("C"));
        assert_eq!(element_symbol(26), Some("Fe"));
        assert_eq!(element_symbol(92), Some("U"));
        assert_eq!(element_symbol(118), Some("Og"));
        assert_eq!(element_symbol(0), None);
        assert_eq!(element_symbol(119), None);
    }

    #[test]
    fn test_atomic_number_from_symbol() {
        assert_eq!(atomic_number_from_symbol("H"), Some(1));
        assert_eq!(atomic_number_from_symbol("h"), Some(1));
        assert_eq!(atomic_number_from_symbol("Fe"), Some(26));
        assert_eq!(atomic_number_from_symbol("fe"), Some(26));
        assert_eq!(atomic_number_from_symbol("FE"), Some(26));
        assert_eq!(atomic_number_from_symbol("U"), Some(92));
        assert_eq!(atomic_number_from_symbol("Og"), Some(118));
        assert_eq!(atomic_number_from_symbol("Xx"), None);
        assert_eq!(atomic_number_from_symbol(""), None);
    }

    #[test]
    fn test_atomic_weight() {
        assert!(atomic_weight(1).unwrap() > 1.0 && atomic_weight(1).unwrap() < 1.1);
        assert!(atomic_weight(6).unwrap() > 12.0 && atomic_weight(6).unwrap() < 12.1);
        assert!(atomic_weight(26).unwrap() > 55.0 && atomic_weight(26).unwrap() < 56.0);
        assert!(atomic_weight(29).unwrap() > 63.0 && atomic_weight(29).unwrap() < 64.0);
        assert_eq!(atomic_weight(0), None);
        assert_eq!(atomic_weight(119), None);
    }

    #[test]
    fn test_k_edge_energy() {
        // Hydrogen and Helium don't have practical K-edges
        assert_eq!(k_edge_energy(1), None);
        assert_eq!(k_edge_energy(2), None);

        // Carbon K-edge should be around 284 eV
        assert!(k_edge_energy(6).unwrap() > 280.0 && k_edge_energy(6).unwrap() < 290.0);

        // Copper K-edge should be around 8979 eV
        assert!(k_edge_energy(29).unwrap() > 8970.0 && k_edge_energy(29).unwrap() < 8990.0);
    }
}
