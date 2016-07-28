import os

chs = ['el','mu']
cats = ['WW','WZ']
binlo	= 900
binhi	= 5000

for ch in chs:

	output = 'docuplots'
	#output = 'docuplots_amATNLO'

	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

	path1 = 'plots_%s_HPW_%s_%s/m_j_fitting/'%(ch,binlo,binhi)
	path2 = 'plots_%s_HPW_%s_%s/m_lvj_fitting/'%(ch,binlo,binhi)

	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< m_j-plots for %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	nameend1 = '%s_mj.pdf'%ch
	print '_TTbar_xwwtreeEDBR_TTBARpowheg_xww_2Gaus_ErfExp_with_pull.pdf'
	os.system('cp %s/_TTbar_xwwtreeEDBR_TTBARpowheg_xww_2Gaus_ErfExp_with_pull.pdf %s/TTbar_%s'%(path1,output,nameend1))
	print '_STop_xwwtreeEDBR_SingleTop_xww_ExpGaus_with_pull.pdf'
	os.system('cp %s/_STop_xwwtreeEDBR_SingleTop_xww_ExpGaus_with_pull.pdf %s/STop_%s'%(path1,output,nameend1))
	print '_VV_xwwtreeEDBR_VV_xww_2_2Gaus_with_pull.pdf'
	os.system('cp %s/_VV_xwwtreeEDBR_VV_xww_2_2Gaus_with_pull.pdf %s/VV_%s'%(path1,output,nameend1))
	print '_WJets0_xwwtreeEDBR_WJets_xww_ErfExp_with_pull.pdf'
	os.system('cp %s/_WJets0_xwwtreeEDBR_WJets_xww_ErfExp_with_pull.pdf %s/WJets_%s'%(path1,output,nameend1))
	print 'm_j_sideband_WJets0_xww__with_pull.pdf'
	os.system('cp %s/m_j_sideband_WJets0_xww__with_pull.pdf %s/bkg_%s'%(path1,output,nameend1))
	print 'm_j_sideband_WJets01_xww__with_pull.pdf'
	os.system('cp %s/m_j_sideband_WJets01_xww__with_pull.pdf %s/bkg_alt_%s'%(path1,output,nameend1))
	print 'm_j_prefit'
	os.system('cp plots_%s_HPW_%s_%s/ExtraPlots/m_j_prefit__WJets0_xww__with_pull.pdf %s/mj_prefit_%s.pdf'%(ch,binlo,binhi,output,ch))

	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< m_lvj-plots for %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	nameend2 = '%s_mlvj_sb.pdf'%ch
	print 'treeEDBR_SingleTop_xww_m_lvj_sb_loExp_with_pull_log.pdf'
	os.system('cp %s/treeEDBR_SingleTop_xww_m_lvj_sb_loExp_with_pull_log.pdf %s/STop_%s'%(path2,output,nameend2))
	print 'treeEDBR_WJets_xww_m_lvj_sb_loExpN_with_pull_log.pdf'
	os.system('cp %s/treeEDBR_WJets_xww_m_lvj_sb_loExpN_with_pull_log.pdf %s/WJets_%s'%(path2,output,nameend2)) 
	print 'treeEDBR_VV_xww_m_lvj_sb_loExp_with_pull_log.pdf'
	os.system('cp %s/treeEDBR_VV_xww_m_lvj_sb_loExp_with_pull_log.pdf %s/VV_%s'%(path2,output,nameend2)) 
	print 'treeEDBR_TTBARpowheg_xww_m_lvj_sb_loExpN_with_pull_log.pdf'
	os.system('cp plots_%s_HPW_%s_%s/ExtraPlots/treeEDBR_TTBARpowheg_xww_m_lvj_sb_loExpN_with_pull_log.pdf %s/TTbar_%s'%(ch,binlo,binhi,output,nameend2))
	print 'm_lvj_sb_lo_WJets0_xww__with_pull_log.pdf'
	os.system('cp plots_%s_HPW_%s_%s/ExtraPlots/m_lvj_sb_lo_WJets0_xww__with_pull_log.pdf %s/all_bkg_%s_mlvj_sb.pdf'%(ch,binlo,binhi,output,ch))




	#print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< closure-test-plots for %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%ch
	#nameend4 = 'closure_%s.pdf'%ch
	#print 'm_j_sideband_WJets0_xww__with_pull.pdf'
	#os.system('cp plots_%s_closure_%s_%s/m_j_fitting/m_j_sideband_WJets0_xww__with_pull.pdf %s/mw_%s'%(ch,binlo,binhi,output,nameend4))
	#print 'm_lvj_sb_lo_WJets0_xww__with_pull_log.pdf'
	#os.system('cp plots_%s_closure_%s_%s/ExtraPlots/m_lvj_sb_lo_WJets0_xww__with_pull_log.pdf %s/sb_%s'%(ch,binlo,binhi,output,nameend4))
	#print 'check_workspace_for_limit__with_pull_log.pdf'
	#os.system('cp plots_%s_closure_%s_%s/m_lvj_fitting/check_workspace_for_limit__with_pull_log.pdf %s/sig_%s'%(ch,binlo,binhi,output,nameend4))


	for cat in cats:

		path3 = 'plots_%s_HP%s_%s_%s/m_lvj_fitting/'%(ch,cat[1],binlo,binhi)
		
		print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< m_lvj-plots in signal region for %s,%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%(ch,cat)
		nameend3 = '%s_mlvj_sig_%s.pdf'%(ch,cat)
		print 'treeEDBR_SingleTop_xww_m_lvj_signal_regionExp_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_SingleTop_xww_m_lvj_signal_regionExp_with_pull_log.pdf %s//STop_%s'%(path3,output,nameend3))
		print 'treeEDBR_WJets_xww_m_lvj_signal_regionExpN_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_WJets_xww_m_lvj_signal_regionExpN_with_pull_log.pdf %s/WJets_%s'%(path3,output,nameend3)) 
		print 'treeEDBR_VV_xww_m_lvj_signal_regionExpN_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_VV_xww_m_lvj_signal_regionExpN_with_pull_log.pdf %s/VV_%s'%(path3,output,nameend3)) 
		print 'treeEDBR_TTBARpowheg_xww_m_lvj_signal_regionExpN_with_pull_log.pdf'
		os.system('cp %s/treeEDBR_TTBARpowheg_xww_m_lvj_signal_regionExpN_with_pull_log.pdf %s/TTbar_%s'%(path3,output,nameend3))

		print 'check_workspace_for_limit__log.pdf'
		os.system('cp plots_%s_HP%s_%s_%s/ExtraPlots/check_workspace_for_limit_log.pdf %s/all_bkg_%s_mlvj_%s.pdf'%(ch,cat[1],binlo,binhi,output,ch,cat))


		print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< alpha-plot for %s,%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'%(ch,cat)
		print 'correction_pdf_WJets0_xww_ExpN_M_lvj_signal_region_to_sideband.pdf'
		os.system('cp plots_%s_HP%s_%s_%s/other/correction_pdf_WJets0_xww_ExpN_M_lvj_signal_region_to_sideband.pdf %s/alpha_%s_%s.pdf'%(ch,cat[1],binlo,binhi,output,ch,cat))

