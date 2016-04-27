import os

chs = ['el','mu']
cats = ['WW','WZ']

for ch in chs:

	#output = 'docuplots'
	output = 'docuplots_amATNLO'

	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

	path1 = 'plots_%s_HPW/m_j_fitting/'%ch
	path2 = 'plots_%s_HPW/m_lvj_fitting/'%ch

	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< m_j-plots for %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	nameend1 = '%s_mj.pdf'%ch
	print '_TTbar_xwwtreeEDBR_TTBARpowheg_xww_2Gaus_ErfExp_with_pull.pdf'
	os.system('cp %s/_TTbar_xwwtreeEDBR_TTBARpowheg_xww_2Gaus_ErfExp_with_pull.pdf %s/mj/TTbar_%s'%(path1,output,nameend1))
	print '_STop_xwwtreeEDBR_SingleTop_xww_ExpGaus_with_pull.pdf'
	os.system('cp %s/_STop_xwwtreeEDBR_SingleTop_xww_ExpGaus_with_pull.pdf %s/mj/STop_%s'%(path1,output,nameend1))
	print '_VV_xwwtreeEDBR_VV_xww_2_2Gaus_with_pull.pdf'
	os.system('cp %s/_VV_xwwtreeEDBR_VV_xww_2_2Gaus_with_pull.pdf %s/mj/VV_%s'%(path1,output,nameend1))
	print '_WJets0_xwwtreeEDBR_WJets_xww_ErfExp_with_pull.pdf'
	os.system('cp %s/_WJets0_xwwtreeEDBR_WJets_xww_ErfExp_with_pull.pdf %s/mj/WJets_%s'%(path1,output,nameend1))
	print 'm_j_sideband_WJets0_xww__with_pull.pdf'
	os.system('cp %s/m_j_sideband_WJets0_xww__with_pull.pdf %s/mj/bkg_%s'%(path1,output,nameend1))

	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< m_lvj-plots for %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	nameend2 = '%s_mlvj_sb.pdf'%ch
	print 'treeEDBR_SingleTop_xww_m_lvj_sb_loExp_with_pull_log.pdf'
	os.system('cp %s/treeEDBR_SingleTop_xww_m_lvj_sb_loExp_with_pull_log.pdf %s/mlvj_sb/STop_%s'%(path2,output,nameend2))
	print 'treeEDBR_WJets_xww_m_lvj_sb_loExpN_with_pull_log.pdf'
	os.system('cp %s/treeEDBR_WJets_xww_m_lvj_sb_loExpN_with_pull_log.pdf %s/mlvj_sb/WJets_%s'%(path2,output,nameend2)) 
	print 'treeEDBR_VV_xww_m_lvj_sb_loExp_with_pull_log.pdf'
	os.system('cp %s/treeEDBR_VV_xww_m_lvj_sb_loExp_with_pull_log.pdf %s/mlvj_sb/VV_%s'%(path2,output,nameend2)) 
	print 'treeEDBR_TTBARpowheg_xww_m_lvj_sb_loExpN_with_pull_log.pdf'
	os.system('cp plots_%s_HPW/ExtraPlots/treeEDBR_TTBARpowheg_xww_m_lvj_sb_loExpN_with_pull_log.pdf %s/mlvj_sb/TTbar_%s'%(ch,output,nameend2))
	print 'm_lvj_sb_lo_WJets0_xww__with_pull_log.pdf'
	os.system('cp plots_%s_HPW/ExtraPlots/m_lvj_sb_lo_WJets0_xww__with_pull_log.pdf %s/mlvj_sb/all_bkg_%s_mlvj_sb.pdf'%(ch,output,ch))




	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< closure-test-plots for %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	nameend4 = 'closure_%s.pdf'%ch
	print 'm_j_sideband_WJets0_xww__with_pull.pdf'
	os.system('cp plots_closure_%s/m_j_fitting/BulkG_WW_lvjj_M2000/m_j_sideband_WJets0_xww__with_pull.pdf %s/closure/mw_%s'%(ch,output,nameend4))
	print 'm_lvj_sb_lo_WJets0_xww__with_pull_log.pdf'
	os.system('cp plots_closure_%s/ExtraPlots/BulkG_WW_lvjj_M2000/m_lvj_sb_lo_WJets0_xww__with_pull_log.pdf %s/closure/sb_%s'%(ch,output,nameend4))
	print 'check_workspace_for_limit__with_pull_log.pdf'
	os.system('cp plots_closure_%s/m_lvj_fitting/BulkG_WW_lvjj_M2000/check_workspace_for_limit__with_pull_log.pdf %s/closure/sig_%s'%(ch,output,nameend4))


	for cat in cats:

		path3 = 'plots_%s_HP%s/m_lvj_fitting/'%(ch,cat[1])
		
		print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< m_lvj-plots in signal region for %s,%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%(ch,cat)
		nameend3 = '%s_mlvj_sig_%s.pdf'%(ch,cat)
		print 'treeEDBR_SingleTop_xww_m_lvj_signal_regionExp_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_SingleTop_xww_m_lvj_signal_regionExp_with_pull_log.pdf %s/mlvj_sig/STop_%s'%(path3,output,nameend3))
		print 'treeEDBR_WJets_xww_m_lvj_signal_regionExpN_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_WJets_xww_m_lvj_signal_regionExpN_with_pull_log.pdf %s/mlvj_sig/WJets_%s'%(path3,output,nameend3)) 
		print 'treeEDBR_VV_xww_m_lvj_signal_regionExpN_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_VV_xww_m_lvj_signal_regionExpN_with_pull_log.pdf %s/mlvj_sig/VV_%s'%(path3,output,nameend3)) 
		print 'treeEDBR_TTBARpowheg_xww_m_lvj_signal_regionExpN_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_TTBARpowheg_xww_m_lvj_signal_regionExpN_with_pull_log.pdf %s/mlvj_sig/TTbar_%s'%(path3,output,nameend3))

		print 'check_workspace_for_limit__log.pdf'
		os.system('cp plots_%s_HP%s/ExtraPlots/check_workspace_for_limit_log.pdf %s/mlvj_sig/all_bkg_%s_mlvj_%s.pdf'%(ch,cat[1],output,ch,cat))


		print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< alpha-plot for %s,%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%(ch,cat)
		print 'correction_pdf_WJets0_xww_ExpN_M_lvj_signal_region_to_sideband.pdf'
		os.system('cp plots_%s_HP%s/other/correction_pdf_WJets0_xww_ExpN_M_lvj_signal_region_to_sideband.pdf %s/mlvj_sig/alpha_%s_%s.pdf'%(ch,cat[1],output,ch,cat))

