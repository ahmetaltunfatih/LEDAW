A. PREPARATION OF SUPERSYSTEM OUTPUT FILE FOR QM/MM INTERACTION ENERGY ANALYSIS

(1) Submit supersystem calculation in the "1ADDUCT" directory.

(2) Submit the calculation in the "2ADDUCT-wo-QM-MM-ELS" reading 
    the "qro" file obtained after running the job in the "1ADDUCT" directory.
    
(3) Substract the printed energy in (2) from the E(0) energy in (1) or from
    the reference energy in the LED part.
    
(4) Insert this energy (QM/MM electrostatics) to the end of the output file in 
    the "1ADDUCT" directory as if it were CPCM dielectric term.
    
    NOTE: This inclusion has been already done. This modified file will serve 
    as the supersystem output file in the LEDAW run.



B. PREPARATION OF SUBSYSTEM OUTPUT FILES FOR QM/MM INTERACTION ENERGY ANALYSIS

(1) Submit subsystem calculations in the "3MONOMERS" directory.

(2) Submit the calculations in the "4MONOMERS-wo-QM-MM-ELS" reading the corresponding
    "qro" files obtained after running the jobs in the "3MONOMERS" directory.
    
(3) Substract the printed energy in (2) from the E(0) energy in (1) for each  
    subsystem.
    
(4) Insert this energy (QM/MM electrostatics) to the end of the corresponding 
    output file in the "3MONOMERS" directory as if it were CPCM dielectric term.
    
    NOTE: This inclusion has been already done for each subsystem. These modified 
    files will serve as the subsystem output files in the LEDAW run.


