import os

chs = ['el','mu']
cats = ['WW','WZ']

for ch in chs:

	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

	path1 = 'plots_%s_HPW/m_j_fitting/BulkG_WW_lvjj_M2000'%ch
	path2 = 'plots_%s_HPW/m_lvj_fitting/BulkG_WW_lvjj_M2000'%ch

	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< m_j-plots for %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	nameend1 = '%s_mj.pdf'%ch
	print '_TTbar_xwwtreeEDBR_TTBARpowheg_xww_2Gaus_ErfExp_with_pull.pdf'
	os.system('cp %s/_TTbar_xwwtreeEDBR_TTBARpowheg_xww_2Gaus_ErfExp_with_pull.pdf docuplots/mj/TTbar_%s'%(path1,nameend1))
	print '_STop_xwwtreeEDBR_SingleTop_xww_ExpGaus_with_pull.pdf'
	os.system('cp %s/_STop_xwwtreeEDBR_SingleTop_xww_ExpGaus_with_pull.pdf docuplots/mj/STop_%s'%(path1,nameend1))
	print '_VV_xwwtreeEDBR_VV_xww_2_2Gaus_with_pull.pdf'
	os.system('cp %s/_VV_xwwtreeEDBR_VV_xww_2_2Gaus_with_pull.pdf docuplots/mj/VV_%s'%(path1,nameend1))
	print '_WJets0_xwwtreeEDBR_WJets_xww_ErfExp_with_pull.pdf'
	os.system('cp %s/_WJets0_xwwtreeEDBR_WJets_xww_ErfExp_with_pull.pdf docuplots/mj/WJets_%s'%(path1,nameend1))
	print 'm_j_sideband_WJets0_xww__with_pull.pdf'
	os.system('cp %s/m_j_sideband_WJets0_xww__with_pull.pdf docuplots/mj/bkg_%s'%(path1,nameend1))

	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< m_lvj-plots for %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	nameend2 = '%s_mlvj_sb.pdf'%ch
	print 'treeEDBR_SingleTop_xww_m_lvj_sb_loExp_with_pull_log.pdf'
	os.system('cp %s/treeEDBR_SingleTop_xww_m_lvj_sb_loExp_with_pull_log.pdf docuplots/mlvj_sb/STop_%s'%(path2,nameend2))
	print 'treeEDBR_WJets_xww_m_lvj_sb_loExpN_with_pull_log.pdf'
	os.system('cp %s/treeEDBR_WJets_xww_m_lvj_sb_loExpN_with_pull_log.pdf docuplots/mlvj_sb/WJets_%s'%(path2,nameend2)) 
	print 'treeEDBR_VV_xww_m_lvj_sb_loExp_with_pull_log.pdf'
	os.system('cp %s/treeEDBR_VV_xww_m_lvj_sb_loExp_with_pull_log.pdf docuplots/mlvj_sb/VV_%s'%(path2,nameend2)) 
	print 'treeEDBR_TTBARpowheg_xww_m_lvj_sb_loExpN_with_pull_log.pdf'
	os.system('cp plots_%s_HPW/ExtraPlots/BulkG_WW_lvjj_M2000/treeEDBR_TTBARpowheg_xww_m_lvj_sb_loExpN_with_pull_log.pdf docuplots/mlvj_sb/TTbar_%s'%(ch,nameend2))
	print 'm_lvj_sb_lo_WJets0_xww__with_pull_log.pdf'
	os.system('cp plots_%s_HPW/ExtraPlots/BulkG_WW_lvjj_M2000/m_lvj_sb_lo_WJets0_xww__with_pull_log.pdf docuplots/mlvj_sb/all_bkg_%s_mlvj_sb.pdf'%(ch,ch))

	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< closure-test-plots for %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	nameend4 = 'closure_%s.pdf'%ch
	print 'm_j_sideband_WJets0_xww__with_pull.pdf'
	os.system('cp plots_closure_%s/m_j_fitting/BulkG_WW_lvjj_M2000/m_j_sideband_WJets0_xww__with_pull.pdf docuplots/closure/mw_%s'%(ch,nameend4))
	print 'm_lvj_sb_lo_WJets0_xww__with_pull_log.pdf'
	os.system('cp plots_closure_%s/ExtraPlots/BulkG_WW_lvjj_M2000/m_lvj_sb_lo_WJets0_xww__with_pull_log.pdf docuplots/closure/sb_%s'%(ch,nameend4))
	print 'check_workspace_for_limit__with_pull_log.pdf'
	os.system('cp plots_closure_%s/m_lvj_fitting/BulkG_WW_lvjj_M2000/check_workspace_for_limit__with_pull_log.pdf docuplots/closure/sig_%s'%(ch,nameend4))


	for cat in cats:

		path3 = 'plots_%s_HP%s/m_lvj_fitting/BulkG_WW_lvjj_M2000'%(ch,cat[1])
		
		print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< m_lvj-plots in signal region for %s,%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%(ch,cat)
		nameend3 = '%s_mlvj_sig_%s.pdf'%(ch,cat)
		print 'treeEDBR_SingleTop_xww_m_lvj_signal_regionExp_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_SingleTop_xww_m_lvj_signal_regionExp_with_pull_log.pdf docuplots/mlvj_sig/STop_%s'%(path3,nameend3))
		print 'treeEDBR_WJets_xww_m_lvj_signal_regionExpN_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_WJets_xww_m_lvj_signal_regionExpN_with_pull_log.pdf docuplots/mlvj_sig/WJets_%s'%(path3,nameend3)) 
		print 'treeEDBR_VV_xww_m_lvj_signal_regionExpN_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_VV_xww_m_lvj_signal_regionExpN_with_pull_log.pdf docuplots/mlvj_sig/VV_%s'%(path3,nameend3)) 
		print 'treeEDBR_TTBARpowheg_xww_m_lvj_signal_regionExpN_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_TTBARpowheg_xww_m_lvj_signal_regionExpN_with_pull_log.pdf docuplots/mlvj_sig/TTbar_%s'%(path3,nameend3))

		print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< alpha-plot for %s,%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%(ch,cat)
		print 'correction_pdf_WJets0_xww_ExpN_M_lvj_signal_region_to_sideband.pdf'
		os.system('cp plots_%s_HP%s/other/BulkG_WW_lvjj_M2000/correction_pdf_WJets0_xww_ExpN_M_lvj_signal_region_to_sideband.pdf docuplots/mlvj_sig/alpha_%s_%s.pdf'%(ch,cat[1],ch,cat))

