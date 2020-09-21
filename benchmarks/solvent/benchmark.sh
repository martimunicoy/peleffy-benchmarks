 
 #PELE with OFF parameters, the OBC implicit solvent and gasteiger charges.
 python main.py OFF -solvent OBC -c gasteiger -p -o OFF_out1

 #PELE with OFF parameters, the OBC implicit solvent and am1bcc charges.
 python main.py OFF -solvent OBC -c am1bcc -p -o OFF_out2

 #PELE with OFF parameters, the OBC implicit solvent and OPLS charges.
 python main.py OFF -solvent OBC -c OPLS -p -o OFF_out3

 #Standard PELE with OPLS2005 parameters and the VDGBNP implicit solvent (use PlopRotTemp in this case)
 python main.py OPLS -solvent VDGNP -p -o OPLS_out4

 #Standard PELE with OPLS2005 parameters and the OBC implicit solvent (use PlopRotTemp in this case)
 python main.py OPLS -solvent OBC -p -o OPLS_out5

 #PELE with OFF bonding parameters, the OPLS2005 non-bonding parameters, VDGBNP implicit solvent and OPLS charges.
 python main.py OFFOPLS -solvent VDGBNP -c OPLS -nb -p -o OFFOPLS_out6

 #PELE with OFF dihedral parameters, the OPLS2005 non-bonding, bond and angle parameters, VDGBNP implicit solvent and OPLS charges.
 python main.py OFFOPLS -solvent VDGBNP -c OPLS -nb -ba -p -o OFFOPLS_out7

 #PELE with OFF bonding parameters, the OPLS2005 non-bonding parameters, VDGBNP implicit solvent and am1bcc charges.
 python main.py OFFOPLS -solvent VDGBNP -c am1bcc -nb -p -o OFFOPLS_out8

 #PELE with OFF dihedral parameters, the OPLS2005 non-bonding, bond and angle parameters, VDGBNP implicit solvent and am1bcc charges.
 python main.py OFFOPLS -solvent VDGBNP -c am1bcc -nb -ba  -p -o OFFOPLS_out9