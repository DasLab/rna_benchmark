#/usr/bin/python

titles = {
	"5P_j12_leadzyme": "5' J1/2, Leadzyme",
	"5P_p1_m_box_riboswitch": "5' P1, M-Box Riboswitch",
	"3P_j55a_group_I_intron": "3' J5/5a, Group I Intron",
	"5P_j55a_group_I_intron": "5' J5/5a, Group I Intron",
	"hepatitis_C_virus_ires_IIa": "Hepatitis C Virus IRES IIa",
	"j24_tpp_riboswitch": "J2/4, TPP Riboswitch",
	"j31_glycine_riboswitch": "J3/1, Glycine Riboswitch",
	"j23_group_II_intron": "J2/3, Group II Intron",
	"l1_sam_II_riboswitch": "L1, SAM-II Riboswitch",
	"l2_viral_rna_pseudoknot": "L2, Viral RNA Pseudoknot",
	"23s_rrna_44_49": "23S rRNA (44-49)",
	"23s_rrna_531_536": "23S rRNA (531-536)",
	"23s_rrna_2534_2540": "23S rRNA (2534-2540)",
	"23s_rrna_1976_1985": "23S rRNA (1976-1985)",
	"23s_rrna_2003_2012": "23S rRNA (2003-2012)"
}

def get_title( target ):
	return titles[target] if target in titles.keys() else target