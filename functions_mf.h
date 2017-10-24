void TMDB::offline_loop( bool ignoreAside )
{
  // offline muon loop
  for(Int_t muon=0; muon<mu_n; muon++) {
    Float_t muon_eta = mu_eta->at(muon);
    Float_t muon_phi = mu_phi->at(muon);
    Float_t muon_pt  = mu_pt->at(muon);
    bool isCombined = (mu_muonType->at(muon)==0)? true: false;
    // bool isCaloTag = (mu_muonType->at(muon)==3)? true: false;
    bool isAside    = (muon_eta>=0)? true: false;
    Int_t quality	  = mu_quality->at(muon);
    // if(muon_pt<15000) {continue;} // pt cut

    th1_offline[0]->Fill(muon_eta);
    if (quality == 0){th1_offline[1]->Fill(muon_eta);} // Tight Working Point
    if (quality == 1){th1_offline[2]->Fill(muon_eta);} // Medium Working Point
    if (quality == 2){th1_offline[3]->Fill(muon_eta);} // Loose Working Point

    // if(!isCombined)   {continue;} // combined muon
    // if(!isCaloTag)   {continue;} // colorimeter tagged  muon

    // Only Combined Muons
    th1_offline_eta[0]->Fill(muon_eta);
    th1_offline_pt[0]->Fill(muon_pt);
    th1_offline_phi[0]->Fill(muon_phi);
    if(isAside){ th1_offline_phi[1]->Fill(muon_phi);}
    else{        th1_offline_phi[5]->Fill(muon_phi);}

    // Working point criterias
    if (quality == 0){ // Tight Working Point
      th1_offline_eta[1]->Fill(muon_eta);
      th1_offline_pt[1]->Fill(muon_pt);
      if(isAside){  th1_offline_phi[2]->Fill(muon_phi);} // A side
      else{         th1_offline_phi[6]->Fill(muon_phi);} // C side
    }
    if (quality == 1){ // Medium Working Point
      th1_offline_eta[2]->Fill(muon_eta);
      th1_offline_pt[2]->Fill(muon_pt);
      if(isAside){  th1_offline_phi[3]->Fill(muon_phi);} // A side
      else{         th1_offline_phi[7]->Fill(muon_phi);} // C side
    }
    if (quality == 2){ // Loose Working Point
      th1_offline_eta[3]->Fill(muon_eta);
      th1_offline_pt[3]->Fill(muon_pt);
      if(isAside){  th1_offline_phi[4]->Fill(muon_phi);} // A side
      else{         th1_offline_phi[8]->Fill(muon_phi);} // C side
    }
  } // end of offline muon loop
} // end of offline loop function

void TMDB::hlt_loop( bool ignoreAside )
{
  // HLT loop
  for(Int_t hlt_loop=0; hlt_loop<trigger_info_n; hlt_loop++) {

    string  chain = trigger_info_chain->at(hlt_loop);
    bool isPassed = (trigger_info_isPassed->at(hlt_loop)==1)? true: false;
    Int_t nTracks = trigger_info_nTracks->at(hlt_loop);

    // cout << "Chain: " << chain << endl;
    if (chain != "HLT_mu26_ivarmedium" && chain != "HLT_mu60" ){continue;} // chain cut
    // if (chain != "HLT_mu26_ivarmedium"){continue;} // chain cut

    // Tracks loop
    // Int_t typeTrack;
    Float_t etaTrack, phiTrack, ptTrack;
    for(Int_t track=0; track<nTracks; track++) { // nTracks loop
      // typeTrack   = trigger_info_typeVec->at(hlt_loop).at(track);
      etaTrack    = trigger_info_etaVec->at(hlt_loop).at(track);
      phiTrack    = trigger_info_phiVec->at(hlt_loop).at(track);
      ptTrack     = trigger_info_ptVec->at(hlt_loop).at(track);

      // if(typeTrack != 0){continue;} // combined muon cut
      bool isAside  = (etaTrack>=0)? true: false;
      bool isBarrel = (fabs(etaTrack)<1.05) ? true: false;

      th1_hlt_eta[0]->Fill(etaTrack);  // Fill HLT eta histogram
      th1_hlt_phi[0]->Fill(phiTrack);  // Fill HLT phi histogram
      if(isAside){  th1_hlt_phi[2]->Fill(phiTrack);} // A side
      else{         th1_hlt_phi[4]->Fill(phiTrack);} // C side
      th1_hlt_pt[0]->Fill(ptTrack);    // Fill HLT pt histogram
      if(isBarrel){ th1_hlt_pt[2]->Fill(ptTrack);} // Barrel region
      else{         th1_hlt_pt[4]->Fill(ptTrack);} // Endcap region

      if(isPassed){ // Fill passed HLT histograms
        th1_hlt_eta[1]->Fill(etaTrack);  // Fill passed HLT eta histogram
        th1_hlt_phi[1]->Fill(phiTrack);  // Fill passed HLT phi histogram
        if(isAside){  th1_hlt_phi[3]->Fill(phiTrack);} // A side
        else{         th1_hlt_phi[5]->Fill(phiTrack);} // C side
        th1_hlt_pt[1]->Fill(ptTrack);    // Fill passed HLT pt histogram
        if(isBarrel){ th1_hlt_pt[3]->Fill(ptTrack);} // Barrel region
        else{         th1_hlt_pt[5]->Fill(ptTrack);} // Endcap region
      }

      // offline muon associated
      for(Int_t muon=0; muon<mu_n; muon++) {
        Float_t muon_eta = mu_eta->at(muon);
        Float_t muon_phi = mu_phi->at(muon);
        Float_t muon_pt  = mu_pt->at(muon);
        Float_t dR = calc_dR(etaTrack-muon_eta , phiTrack-muon_phi);
        // bool isCombined = (mu_muonType->at(muon)==0)? true: false;
        // bool isCaloTag = (mu_muonType->at(muon)==3)? true: false;
        Int_t quality	  = mu_quality->at(muon);
        // if(!isCombined)   {continue;} // combined muon cut
        // if(dR> 0.1)       {continue;} // dR cut

        // true positive: it has an offline muon associated
        if(isPassed && dR < 0.1){
          th1_hlt_eta_vp[0]->Fill(etaTrack);
          th1_hlt_phi_vp[0]->Fill(phiTrack);
          if(isAside){  th1_hlt_phi_vp[1]->Fill(phiTrack);} // A side
          else{         th1_hlt_phi_vp[5]->Fill(phiTrack);} // C side
          th1_hlt_pt_vp[0]->Fill(ptTrack);
          if(isBarrel){ th1_hlt_pt_vp[4]->Fill(ptTrack);} // Barrel region
          else{         th1_hlt_pt_vp[8]->Fill(ptTrack);} // Endcap region
          if(quality==0){ // Tight Working Point
            th1_hlt_eta_vp[1]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_vp[2]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_vp[6]->Fill(phiTrack);} // C side
            th1_hlt_pt_vp[1]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_vp[5]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_vp[9]->Fill(ptTrack);} // Endcap region
          }
          if(quality==1){ // Medium Working Point
            th1_hlt_eta_vp[2]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_vp[3]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_vp[7]->Fill(phiTrack);} // C side
            th1_hlt_pt_vp[2]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_vp[6]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_vp[10]->Fill(ptTrack);} // Endcap region
          }
          if(quality==2){ // Loose Working Point
            th1_hlt_eta_vp[3]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_vp[4]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_vp[8]->Fill(phiTrack);} // C side
            th1_hlt_pt_vp[3]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_vp[7]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_vp[11]->Fill(ptTrack);} // Endcap region
          }
          break; // break loop: muon already found!
        } // end of true positive condition

        // false negative: it doesn't have an offline muon associated
        if(!isPassed && dR < 0.1){
          th1_hlt_eta_fn[0]->Fill(etaTrack);
          th1_hlt_phi_fn[0]->Fill(phiTrack);
          if(isAside){  th1_hlt_phi_fn[1]->Fill(phiTrack);} // A side
          else{         th1_hlt_phi_fn[5]->Fill(phiTrack);} // C side
          th1_hlt_pt_fn[0]->Fill(ptTrack);
          if(isBarrel){ th1_hlt_pt_fn[4]->Fill(ptTrack);} // Barrel region
          else{         th1_hlt_pt_fn[8]->Fill(ptTrack);} // Endcap region
          if(quality==0){ // Tight Working Point
            th1_hlt_eta_fn[1]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_fn[2]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_fn[6]->Fill(phiTrack);} // C side
            th1_hlt_pt_fn[1]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_fn[5]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_fn[9]->Fill(ptTrack);} // Endcap region
          }
          if(quality==1){ // Medium Working Point
            th1_hlt_eta_fn[2]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_fn[3]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_fn[7]->Fill(phiTrack);} // C side
            th1_hlt_pt_fn[2]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_fn[6]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_fn[10]->Fill(ptTrack);} // Endcap region
          }
          if(quality==2){ // Loose Working Point
            th1_hlt_eta_fn[3]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_fn[4]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_fn[8]->Fill(phiTrack);} // C side
            th1_hlt_pt_fn[3]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_fn[7]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_fn[11]->Fill(ptTrack);} // Endcap region
          }
          break; // break loop: muon already found!
        } // end of false negative condition

        // false positive: it has an offline muon associated
        if(isPassed && dR > 0.1 && muon==mu_n-1){
          th1_hlt_eta_fp[0]->Fill(etaTrack);
          th1_hlt_phi_fp[0]->Fill(phiTrack);
          if(isAside){  th1_hlt_phi_fp[1]->Fill(phiTrack);} // A side
          else{         th1_hlt_phi_fp[5]->Fill(phiTrack);} // C side
          th1_hlt_pt_fp[0]->Fill(ptTrack);
          if(isBarrel){ th1_hlt_pt_fp[4]->Fill(ptTrack);} // Barrel region
          else{         th1_hlt_pt_fp[8]->Fill(ptTrack);} // Endcap region
          if(quality==0){ // Tight Working Point
            th1_hlt_eta_fp[1]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_fp[2]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_fp[6]->Fill(phiTrack);} // C side
            th1_hlt_pt_fp[1]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_fp[5]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_fp[9]->Fill(ptTrack);} // Endcap region
          }
          if(quality==1){ // Medium Working Point
            th1_hlt_eta_fp[2]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_fp[3]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_fp[7]->Fill(phiTrack);} // C side
            th1_hlt_pt_fp[2]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_fp[6]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_fp[10]->Fill(ptTrack);} // Endcap region
          }
          if(quality==2){ // Loose Working Point
            th1_hlt_eta_fp[3]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_fp[4]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_fp[8]->Fill(phiTrack);} // C side
            th1_hlt_pt_fp[3]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_fp[7]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_fp[11]->Fill(ptTrack);} // Endcap region
          }
        } // end of false positive condition

        // true negative: it doesn't have an offline muon associated
        if(!isPassed && dR > 0.1 && muon==mu_n-1){
          th1_hlt_eta_vn[0]->Fill(etaTrack);
          th1_hlt_phi_vn[0]->Fill(phiTrack);
          if(isAside){  th1_hlt_phi_vn[1]->Fill(phiTrack);} // A side
          else{         th1_hlt_phi_vn[5]->Fill(phiTrack);} // C side
          th1_hlt_pt_vn[0]->Fill(ptTrack);
          if(isBarrel){ th1_hlt_pt_vn[4]->Fill(ptTrack);} // Barrel region
          else{         th1_hlt_pt_vn[8]->Fill(ptTrack);} // Endcap region
          if(quality==0){ // Tight Working Point
            th1_hlt_eta_vn[1]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_vn[2]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_vn[6]->Fill(phiTrack);} // C side
            th1_hlt_pt_vn[1]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_vn[5]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_vn[9]->Fill(ptTrack);} // Endcap region
          }
          if(quality==1){ // Medium Working Point
            th1_hlt_eta_vn[2]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_vn[3]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_vn[7]->Fill(phiTrack);} // C side
            th1_hlt_pt_vn[2]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_vn[6]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_vn[10]->Fill(ptTrack);} // Endcap region
          }
          if(quality==2){ // Loose Working Point
            th1_hlt_eta_vn[3]->Fill(etaTrack);
            if(isAside){  th1_hlt_phi_vn[4]->Fill(phiTrack);} // A side
            else{         th1_hlt_phi_vn[8]->Fill(phiTrack);} // C side
            th1_hlt_pt_vn[3]->Fill(ptTrack);
            if(isBarrel){ th1_hlt_pt_vn[7]->Fill(ptTrack);} // Barrel region
            else{         th1_hlt_pt_vn[11]->Fill(ptTrack);} // Endcap region
          }
        } // end of true negative condition
      } // end of muon loop
    } // end of nTracks loop
  } // end of HLT loop
} // end of HLT loop funtion

void TMDB::l1_trigger_loop( bool ignoreAside )
{
  // L1 Trigger loop
  for(Int_t nloop=0; nloop<trig_L1_mu_n; nloop++) {
    Float_t eta   = trig_L1_mu_eta->at(nloop);
    Float_t phi   = trig_L1_mu_phi->at(nloop);
    Int_t pt_thr  = trig_L1_mu_thrNumber->at(nloop);
    Int_t roi     = trig_L1_mu_RoINumber->at(nloop);
    Int_t ssc     = (roi+4)/8;
    Int_t sector  = (trig_L1_mu_sectorAddress->at(nloop)>>1) & (0x3F);
    Int_t vetoed  = trig_L1_mu_vetoed->at(nloop);
    bool isAside  = (trig_L1_mu_hemisphere->at(nloop)==1) ? true: false;
    // bool isAside  = (eta>0.0) ? true: false;
    bool isEndcap = (trig_L1_mu_source->at(nloop)==1) ? true: false;
    bool isBarrel = (trig_L1_mu_source->at(nloop)==0) ? true: false;
    Int_t tile_module = tilemodule_to_check(sector);
    // if(!(isEndcap && pt_thr >= 5 && ssc < 6)) { continue; }
    // if( fabs(eta) < 1.0 || 1.3 < fabs(eta) )  { continue; }
    // if(vetoed) {continue;}
    if( ignoreAside && isAside) {continue;} // test

    // Filling L1 histograms
    th1_l1muon->Fill(eta);
    // pt_thr: 5 = L1MU15-L1MU20, >=5 = L1MU15, 6 = L1MU20
    // For L1MU15-L1MU20
    if (pt_thr == 5) {th1_l1muon_eta->Fill(eta); th1_l1muon_phi->Fill(phi);}
    // For L1MU15
    if (pt_thr >= 5) {
      th1_l1muon_15_eta->Fill(eta);
      th1_l1muon_15_phi[0]->Fill(phi);
      if(isAside){  th1_l1muon_15_phi[1]->Fill(phi);}
      else{         th1_l1muon_15_phi[2]->Fill(phi);}
    }
    // For L1MU20
    if (pt_thr == 6) {
      th1_l1muon_20_eta->Fill(eta);
      th1_l1muon_20_phi[0]->Fill(phi);
      if(isAside){  th1_l1muon_20_phi[1]->Fill(phi);}
      else{         th1_l1muon_20_phi[2]->Fill(phi);}
    }

    // offline muon associated
    for(Int_t muon=0; muon<mu_n; muon++) {
      Float_t muon_eta = mu_eta->at(muon);
      Float_t muon_phi = mu_phi->at(muon);
      Float_t muon_pt  = mu_pt->at(muon);
      // bool isCombined = (mu_muonType->at(muon)==0)? true: false;
      Int_t quality	  = mu_quality->at(muon);

      Float_t dR = calc_dR(eta-muon_eta , phi-muon_phi);

      // true positive: it has an offline muon associated
      if(dR < 0.1){

      } // end of true positive

    } // end of muon loop

    // if( fabs(eta) < 1.0 || 1.3 < fabs(eta) )  { continue; } // eta restriction for TMDB
    // // Get data from TMDB
    // for(Int_t tmdb_n=0; tmdb_n<TILE_murcv_trig_n; tmdb_n++){
    // 	if (!(tile_module == TILE_murcv_trig_mod->at(tmdb_n))) {continue;} // tile module different
    //   if(TILE_murcv_trig_bit0->at(tmdb_n)==1 || TILE_murcv_trig_bit2->at(tmdb_n)==1) {
    //     th1_tmdb->Fill(eta);
    //   }
    // } // end of tmdb coincidence loop
  } // end of L1 Trigger nloop

} // End of tgc_trigger_loop function

// ###############################################################################
// Int_t TMDB::tile_phi_module(Float_t phi)
// Description: Returns the Tile module associated with phi value
// Requirements: none
// ###############################################################################
Int_t TMDB::tile_phi_module(Float_t phi){

	Int_t module = -1;

	Float_t sector	= TMath::Pi()/32.0;
	Int_t	  pos		= fabs(phi)/sector;

	if(phi>=0){ module = pos;}
	else{module = 63 - pos;}

	return module;
} // End of tile_phi_module

// function to print histograms in a ROOT File
void TMDB::print_root(string fout)
{

  string rootfile = fout + ".root";
  TFile *output = new TFile(rootfile.c_str(), "recreate");

  // Writing offline histograms
  for (Int_t i=0; i<4; i++){th1_offline[i]->Write();}
  for (Int_t i=0; i<4; i++){th1_offline_eta[i]->Write();}
  for (Int_t i=0; i<9; i++){th1_offline_phi[i]->Write();}
  for (Int_t i=0; i<4; i++){th1_offline_pt[i]->Write();}

  // Writing HLT histograms
  for (Int_t i=0; i<2; i++){  th1_hlt_eta[i]->Write();}
  for (Int_t i=0; i<4; i++){  th1_hlt_eta_vp[i]->Write();}
  for (Int_t i=0; i<4; i++){  th1_hlt_eta_vn[i]->Write();}
  for (Int_t i=0; i<4; i++){  th1_hlt_eta_fp[i]->Write();}
  for (Int_t i=0; i<4; i++){  th1_hlt_eta_fn[i]->Write();}
  for (Int_t i=0; i<6; i++){  th1_hlt_phi[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_hlt_phi_vp[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_hlt_phi_vn[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_hlt_phi_fp[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_hlt_phi_fn[i]->Write();}
  for (Int_t i=0; i<6; i++){  th1_hlt_pt[i]->Write();}
  for (Int_t i=0; i<12; i++){ th1_hlt_pt_vp[i]->Write();}
  for (Int_t i=0; i<12; i++){ th1_hlt_pt_vn[i]->Write();}
  for (Int_t i=0; i<12; i++){ th1_hlt_pt_fp[i]->Write();}
  for (Int_t i=0; i<12; i++){ th1_hlt_pt_fn[i]->Write();}

  // Writing  L1 histograms
  th1_l1muon->Write();
  th1_l1muon_eta->Write();
  th1_l1muon_phi->Write();
  th1_l1muon_15_eta->Write();
  for (Int_t i=0; i<4; i++){  th1_l1muon_15_eta_vp[i]->Write();}
  for (Int_t i=0; i<4; i++){  th1_l1muon_15_eta_vn[i]->Write();}
  for (Int_t i=0; i<4; i++){  th1_l1muon_15_eta_fp[i]->Write();}
  for (Int_t i=0; i<4; i++){  th1_l1muon_15_eta_fn[i]->Write();}
  for (Int_t i=0; i<3; i++){  th1_l1muon_15_phi[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_15_phi_vp[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_15_phi_vn[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_15_phi_fp[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_15_phi_fn[i]->Write();}
  for (Int_t i=0; i<3; i++){  th1_l1muon_15_pt[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_15_pt_vp[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_15_pt_vn[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_15_pt_fp[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_15_pt_fn[i]->Write();}
  th1_l1muon_20_eta->Write();
  for (Int_t i=0; i<4; i++){  th1_l1muon_20_eta_vp[i]->Write();}
  for (Int_t i=0; i<4; i++){  th1_l1muon_20_eta_vn[i]->Write();}
  for (Int_t i=0; i<4; i++){  th1_l1muon_20_eta_fp[i]->Write();}
  for (Int_t i=0; i<4; i++){  th1_l1muon_20_eta_fn[i]->Write();}
  for (Int_t i=0; i<3; i++){  th1_l1muon_20_phi[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_20_phi_vp[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_20_phi_vn[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_20_phi_fp[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_20_phi_fn[i]->Write();}
  for (Int_t i=0; i<3; i++){  th1_l1muon_20_pt[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_20_pt_vp[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_20_pt_vn[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_20_pt_fp[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_20_pt_fn[i]->Write();}


  // th1_tmdb->Write();

  output->Write();
  output->Close();

}

Float_t TMDB::calc_dR(Float_t dEta, Float_t dPhi)
{
  Float_t dR = sqrt( (dEta*dEta) + (dPhi*dPhi) );
  return(dR);
}

Int_t TMDB::tilemodule_to_check(Int_t sector)
{
  Int_t module = 0;
  switch ( sector )
    {
    case 0:
      module = 61;
      break;
    case 1:
      module = 62;
      break;

    case 2:
      module =  0;
      break;
    case 3:
      module =  1;
      break;
    case 4:
      module =  2;
      break;

    case 5:
      module =  4;
      break;
    case 6:
      module =  5;
      break;
    case 7:
      module =  6;
      break;

    case 8:
      module =  8;
      break;
    case 9:
      module =  9;
      break;
    case 10:
      module = 10;
      break;

    case 11:
      module = 12;
      break;
    case 12:
      module = 13;
      break;
    case 13:
      module = 14;
      break;

    case 14:
      module = 16;
      break;
    case 15:
      module = 17;
      break;
    case 16:
      module = 18;
      break;

    case 17:
      module = 20;
      break;
    case 18:
      module = 21;
      break;
    case 19:
      module = 22;
      break;

    case 20:
      module = 24;
      break;
    case 21:
      module = 25;
      break;
    case 22:
      module = 26;
      break;

    case 23:
      module = 28;
      break;
    case 24:
      module = 29;
      break;
    case 25:
      module = 30;
      break;

    case 26:
      module = 32;
      break;
    case 27:
      module = 33;
      break;
    case 28:
      module = 34;
      break;

    case 29:
      module = 36;
      break;
    case 30:
      module = 37;
      break;
    case 31:
      module = 38;
      break;

    case 32:
      module = 40;
      break;
    case 33:
      module = 41;
      break;
    case 34:
      module = 42;
      break;

    case 35:
      module = 44;
      break;
    case 36:
      module = 45;
      break;
    case 37:
      module = 46;
      break;

    case 38:
      module = 48;
      break;
    case 39:
      module = 49;
      break;
    case 40:
      module = 50;
      break;

    case 41:
      module = 52;
      break;
    case 42:
      module = 53;
      break;
    case 43:
      module = 54;
      break;

    case 44:
      module = 56;
      break;
    case 45:
      module = 57;
      break;
    case 46:
      module = 58;
      break;

    case 47:
      module = 60;
      break;
    }

  return module;
} // End of tilemodule_to_check
