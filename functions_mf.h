// void TMDB::loop_mc( bool ignoreAside )
// {
//
//   // Monte carlo truth loop
//   for(Int_t muon=0; muon<mc_n; muon++) {
//
//   }
// }
void TMDB::loop_offline( bool ignoreAside )
{
  // offline muon loop
  for(Int_t muon=0; muon<mu_n; muon++) {
    Float_t muon_eta = mu_eta->at(muon);
    Float_t muon_phi = mu_phi->at(muon);
    Float_t muon_pt  = mu_pt->at(muon);
    // bool isCombined = (mu_muonType->at(muon)==0)? true: false;
    bool isBarrel = (fabs(muon_eta)<1.05)? true: false;
    bool isEndcap = (fabs(muon_eta)>1.05 && fabs(muon_eta)<2.4)? true: false;
    bool isTMDBroi= (fabs(muon_eta)>0.95 && fabs(muon_eta)<1.3)? true: false;
    bool isCaloTag = (mu_muonType->at(muon)==3)? true: false;
    bool isAside    = (muon_eta>=0)? true: false;
    Int_t quality	  = mu_quality->at(muon);
    if(isCombined) {continue;} // muon combined cut
    // if(muon_pt<10000) {continue;} // pt cut

    if(muon_pt>20000) {th1_offline_eta[4]->Fill(muon_eta);} // pt cut for 20gev
    if(isCaloTag){th1_offline_eta[8]->Fill(muon_eta);} // calorimeter tagged muon
    if(isTMDBroi){th1_offline_pt[12]->Fill(muon_pt);} // TMDB region

    // All muons
    th1_offline_eta[0]->Fill(muon_eta);
    th1_offline_phi[0]->Fill(muon_phi);
    if(isAside){ th1_offline_phi[4]->Fill(muon_phi);}
    else{        th1_offline_phi[8]->Fill(muon_phi);}
    th1_offline_pt[0]->Fill(muon_pt);
    if(isBarrel){th1_offline_pt[4]->Fill(muon_pt);}
    if(isEndcap){th1_offline_pt[8]->Fill(muon_pt);}

    // ------ Working point criterias -------
    // Tight Working Point
    if (quality == 0){
      // Eta Filling
      th1_offline_eta[1]->Fill(muon_eta); // Tight hist
      th1_offline_eta[2]->Fill(muon_eta); // Medium hist
      th1_offline_eta[3]->Fill(muon_eta); // Loose hist
      // Phi Filling
      th1_offline_phi[1]->Fill(muon_phi); // Tight hist
      th1_offline_phi[2]->Fill(muon_phi); // Medium hist
      th1_offline_phi[3]->Fill(muon_phi); // Loose hist
      if(isAside){ // A side
        th1_offline_phi[5]->Fill(muon_phi); // Tight hist
        th1_offline_phi[6]->Fill(muon_phi); // Medium hist
        th1_offline_phi[7]->Fill(muon_phi); // Loose hist
      }else{ // C side
        th1_offline_phi[9]->Fill(muon_phi); // Tight hist
        th1_offline_phi[10]->Fill(muon_phi); // Medium hist
        th1_offline_phi[11]->Fill(muon_phi); // Loose hist
      }
      // Pt Filling
      th1_offline_pt[1]->Fill(muon_pt); // Tight hist
      th1_offline_pt[2]->Fill(muon_pt); // Medium hist
      th1_offline_pt[3]->Fill(muon_pt); // Loose hist
      if(muon_pt>20000) {
        th1_offline_eta[5]->Fill(muon_eta); // pt cut for 20gev
        th1_offline_eta[6]->Fill(muon_eta); // pt cut for 20gev
        th1_offline_eta[7]->Fill(muon_eta); // pt cut for 20gev
      }
      if(isTMDBroi){
        th1_offline_pt[13]->Fill(muon_pt);
        th1_offline_pt[14]->Fill(muon_pt);
        th1_offline_pt[15]->Fill(muon_pt);
      }
      if(isBarrel){ // Barrel region
        th1_offline_pt[5]->Fill(muon_pt); // Tight hist
        th1_offline_pt[6]->Fill(muon_pt); // Medium hist
        th1_offline_pt[7]->Fill(muon_pt); // Loose hist
      }
      if(isEndcap){ // Endcap region
        th1_offline_pt[9]->Fill(muon_pt); // Tight hist
        th1_offline_pt[10]->Fill(muon_pt); // Medium hist
        th1_offline_pt[11]->Fill(muon_pt); // Loose hist
      }
    } // end of Tight Working Point

    // Medium Working Point
    if (quality == 1){
      // Eta Filling
      th1_offline_eta[2]->Fill(muon_eta); // Medium hist
      th1_offline_eta[3]->Fill(muon_eta); // Loose hist
      // Phi Filling
      th1_offline_phi[2]->Fill(muon_phi); // Medium hist
      th1_offline_phi[3]->Fill(muon_phi); // Loose hist
      if(isAside){ // A side
        th1_offline_phi[6]->Fill(muon_phi); // Medium hist
        th1_offline_phi[7]->Fill(muon_phi); // Loose hist
      }else{ // C side
        th1_offline_phi[10]->Fill(muon_phi); // Medium hist
        th1_offline_phi[11]->Fill(muon_phi); // Loose hist
      }
      // Pt Filling
      th1_offline_pt[2]->Fill(muon_pt); // Medium hist
      th1_offline_pt[3]->Fill(muon_pt); // Loose hist
      if(muon_pt>20000) {
        th1_offline_eta[6]->Fill(muon_eta); // pt cut for 20gev
        th1_offline_eta[7]->Fill(muon_eta); // pt cut for 20gev
      }
      if(isTMDBroi){
        th1_offline_pt[14]->Fill(muon_pt);
        th1_offline_pt[15]->Fill(muon_pt);
      }
      if(isBarrel){ // Barrel region
        th1_offline_pt[6]->Fill(muon_pt); // Medium hist
        th1_offline_pt[7]->Fill(muon_pt); // Loose hist
      }
      if(isEndcap){ // Endcap region
        th1_offline_pt[10]->Fill(muon_pt); // Medium hist
        th1_offline_pt[11]->Fill(muon_pt); // Loose hist
      }
    } // enf of Medium Working Point

    // Loose Working Point
    if (quality == 2){
      // Eta Filling
      th1_offline_eta[3]->Fill(muon_eta); // Loose hist
      // Phi Filling
      th1_offline_phi[3]->Fill(muon_phi); // Loose hist
      if(isAside){ // A side
        th1_offline_phi[7]->Fill(muon_phi); // Loose hist
      }else{ // C side
        th1_offline_phi[11]->Fill(muon_phi); // Loose hist
      }
      // Pt Filling
      th1_offline_pt[3]->Fill(muon_pt); // Loose hist
      if(muon_pt>20000) {
        th1_offline_eta[7]->Fill(muon_eta); // pt cut for 20gev
      }
      if(isTMDBroi){th1_offline_pt[15]->Fill(muon_pt);}
      if(isBarrel){ // Barrel region
        th1_offline_pt[7]->Fill(muon_pt); // Loose hist
      }
      if(isEndcap){ // Endcap region
        th1_offline_pt[11]->Fill(muon_pt); // Loose hist
      }
    } // end of Loose Working Point
    // ------ End of Working point criterias -------

  } // end of offline muon loop
} // end of offline loop function

void TMDB::loop_hlt(bool ignoreAside)
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
        // Float_t muon_pt  = mu_pt->at(muon);
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

void TMDB::loop_l1_trigger( bool ignoreAside )
{
  // L1 Trigger loop
  for(Int_t nloop=0; nloop<trig_L1_mu_n; nloop++) {
    Float_t eta   = trig_L1_mu_eta->at(nloop);
    Float_t phi   = trig_L1_mu_phi->at(nloop);
    Int_t pt_thr  = trig_L1_mu_thrNumber->at(nloop);
    // Int_t roi     = trig_L1_mu_RoINumber->at(nloop);
    // Int_t ssc     = (roi+4)/8;
    // Int_t sector  = (trig_L1_mu_sectorAddress->at(nloop)>>1) & (0x3F);
    // Int_t tile_module = tilemodule_to_check(sector);
    // Int_t vetoed  = trig_L1_mu_vetoed->at(nloop);
    bool isAside  = (trig_L1_mu_hemisphere->at(nloop)==1) ? true: false;
    // bool isAside  = (eta>0.0) ? true: false;
    bool isBarrel = (trig_L1_mu_source->at(nloop)==0) ? true: false;
    bool isEndcap = (trig_L1_mu_source->at(nloop)==1) ? true: false;
    // bool isForward= (trig_L1_mu_source->at(nloop)==2) ? true: false;
    // if(!(isEndcap && pt_thr >= 5 && ssc < 6)) { continue; }
    // if( fabs(eta) < 1.0 || 1.3 < fabs(eta) )  { continue; }
    // if(vetoed) {continue;}
    // if( ignoreAside && isAside) {continue;} // test
    // if( isForward) {continue;} // test

    // Filling L1 histograms
    th1_l1muon->Fill(eta);
    // pt_thr: 5 = L1MU15-L1MU20, >=5 = L1MU15, 6 = L1MU20
    // For L1MU15-L1MU20
    if (pt_thr == 5) {th1_l1muon_eta->Fill(eta); th1_l1muon_phi->Fill(phi);}
    // For L1MU15
    if (pt_thr >= 5) {
      th1_l1muon_15_eta->Fill(eta);
      th1_l1muon_15_phi[0]->Fill(phi);
      if(isBarrel){
        th1_l1muon_15_phi[1]->Fill(phi);
        if(isAside){  th1_l1muon_15_phi[4]->Fill(phi);}
        else{         th1_l1muon_15_phi[7]->Fill(phi);}
      }
      if(isEndcap){
        th1_l1muon_15_phi[2]->Fill(phi);
        if(isAside){  th1_l1muon_15_phi[5]->Fill(phi);}
        else{         th1_l1muon_15_phi[8]->Fill(phi);}
      }
      if(isAside){  th1_l1muon_15_phi[3]->Fill(phi);}
      else{         th1_l1muon_15_phi[6]->Fill(phi);}
    } // end of L1MU15 if

    // For L1MU20
    if (pt_thr >= 6) {
      th1_l1muon_20_eta->Fill(eta);
      th1_l1muon_20_phi[0]->Fill(phi);
      if(isBarrel){
        th1_l1muon_20_phi[1]->Fill(phi);
        if(isAside){  th1_l1muon_20_phi[4]->Fill(phi);}
        else{         th1_l1muon_20_phi[7]->Fill(phi);}
      }
      if(isEndcap){
        th1_l1muon_20_phi[2]->Fill(phi);
        if(isAside){  th1_l1muon_20_phi[5]->Fill(phi);}
        else{         th1_l1muon_20_phi[8]->Fill(phi);}
      }
      if(isAside){  th1_l1muon_20_phi[3]->Fill(phi);}
      else{         th1_l1muon_20_phi[6]->Fill(phi);}
    } // end of L1MU20 if

    // offline muon associated
    for(Int_t muon=0; muon<mu_n; muon++) {
      Float_t muon_eta = mu_eta->at(muon);
      Float_t muon_phi = mu_phi->at(muon);
      Float_t muon_pt  = mu_pt->at(muon);
      // bool isCombined = (mu_muonType->at(muon)==0)? true: false;
      Int_t quality	  = mu_quality->at(muon);

      Float_t dR = calc_dR(eta-muon_eta , phi-muon_phi);

      // Total offline muon pt / L1MU15
      if (pt_thr >= 5){
        th1_l1muon_15_pt[0]->Fill(muon_pt);
        if(isBarrel){th1_l1muon_15_pt[1]->Fill(muon_pt);}
        if(isEndcap){th1_l1muon_15_pt[2]->Fill(muon_pt);}
      }
      // Total offline muon pt / L1MU20
      if (pt_thr >= 6){
        th1_l1muon_20_pt[0]->Fill(muon_pt);
        if(isBarrel){th1_l1muon_20_pt[1]->Fill(muon_pt);}
        if(isEndcap){th1_l1muon_20_pt[2]->Fill(muon_pt);}
      }

      // true positive: it has an offline muon associated
      if(dR < 0.1){
        // for L1MU15
        if (pt_thr >= 5){
          th1_l1muon_15_eta_vp[0]->Fill(eta);
          th1_l1muon_15_phi_vp[0]->Fill(phi);
          th1_l1muon_15_pt_vp[0]->Fill(muon_pt);
          if(isBarrel){
            th1_l1muon_15_phi_vp[4]->Fill(phi);
            th1_l1muon_15_pt_vp[4]->Fill(muon_pt);
            if(isAside){  th1_l1muon_15_phi_vp[16]->Fill(phi);}
            else{         th1_l1muon_15_phi_vp[28]->Fill(phi);}
          }
          if(isEndcap){
            th1_l1muon_15_phi_vp[8]->Fill(phi);
            th1_l1muon_15_pt_vp[8]->Fill(muon_pt);
            if(isAside){  th1_l1muon_15_phi_vp[20]->Fill(phi);}
            else{         th1_l1muon_15_phi_vp[32]->Fill(phi);}
          }
          if(isAside){  th1_l1muon_15_phi_vp[12]->Fill(phi);}
          else{         th1_l1muon_15_phi_vp[24]->Fill(phi);}

          // Quality criteria
          switch (quality) {
            case 0 /* tight */ :{
              th1_l1muon_15_eta_vp[1]->Fill(eta);
              th1_l1muon_15_eta_vp[2]->Fill(eta);
              th1_l1muon_15_eta_vp[3]->Fill(eta);
              th1_l1muon_15_phi_vp[1]->Fill(phi);
              th1_l1muon_15_phi_vp[2]->Fill(phi);
              th1_l1muon_15_phi_vp[3]->Fill(phi);
              if(isAside){
                th1_l1muon_15_phi_vp[13]->Fill(phi);
                th1_l1muon_15_phi_vp[14]->Fill(phi);
                th1_l1muon_15_phi_vp[15]->Fill(phi);
              }else{
                th1_l1muon_15_phi_vp[25]->Fill(phi);
                th1_l1muon_15_phi_vp[26]->Fill(phi);
                th1_l1muon_15_phi_vp[27]->Fill(phi);
              }
              th1_l1muon_15_pt_vp[1]->Fill(muon_pt);
              th1_l1muon_15_pt_vp[2]->Fill(muon_pt);
              th1_l1muon_15_pt_vp[3]->Fill(muon_pt);
              if(isBarrel){
                th1_l1muon_15_pt_vp[5]->Fill(muon_pt);
                th1_l1muon_15_pt_vp[6]->Fill(muon_pt);
                th1_l1muon_15_pt_vp[7]->Fill(muon_pt);
                th1_l1muon_15_phi_vp[5]->Fill(phi);
                th1_l1muon_15_phi_vp[6]->Fill(phi);
                th1_l1muon_15_phi_vp[7]->Fill(phi);
                if(isAside){
                  th1_l1muon_15_phi_vp[17]->Fill(phi);
                  th1_l1muon_15_phi_vp[18]->Fill(phi);
                  th1_l1muon_15_phi_vp[19]->Fill(phi);
                }else{
                  th1_l1muon_15_phi_vp[29]->Fill(phi);
                  th1_l1muon_15_phi_vp[30]->Fill(phi);
                  th1_l1muon_15_phi_vp[31]->Fill(phi);
                }
              }
              if(isEndcap){
                th1_l1muon_15_pt_vp[9]->Fill(muon_pt);
                th1_l1muon_15_pt_vp[10]->Fill(muon_pt);
                th1_l1muon_15_pt_vp[11]->Fill(muon_pt);
                th1_l1muon_15_phi_vp[9]->Fill(phi);
                th1_l1muon_15_phi_vp[10]->Fill(phi);
                th1_l1muon_15_phi_vp[11]->Fill(phi);
                if(isAside){
                  th1_l1muon_15_phi_vp[21]->Fill(phi);
                  th1_l1muon_15_phi_vp[22]->Fill(phi);
                  th1_l1muon_15_phi_vp[23]->Fill(phi);
                }
                else{
                  th1_l1muon_15_phi_vp[33]->Fill(phi);
                  th1_l1muon_15_phi_vp[34]->Fill(phi);
                  th1_l1muon_15_phi_vp[35]->Fill(phi);
                }
              }
              break;
            }
            case 1 /* medium */:{
              th1_l1muon_15_eta_vp[2]->Fill(eta);
              th1_l1muon_15_eta_vp[3]->Fill(eta);
              th1_l1muon_15_phi_vp[2]->Fill(phi);
              th1_l1muon_15_phi_vp[3]->Fill(phi);
              if(isAside){
                th1_l1muon_15_phi_vp[14]->Fill(phi);
                th1_l1muon_15_phi_vp[15]->Fill(phi);
              }else{
                th1_l1muon_15_phi_vp[26]->Fill(phi);
                th1_l1muon_15_phi_vp[27]->Fill(phi);
              }
              th1_l1muon_15_pt_vp[2]->Fill(muon_pt);
              th1_l1muon_15_pt_vp[3]->Fill(muon_pt);
              if(isBarrel){
                th1_l1muon_15_pt_vp[6]->Fill(muon_pt);
                th1_l1muon_15_pt_vp[7]->Fill(muon_pt);
                th1_l1muon_15_phi_vp[6]->Fill(phi);
                th1_l1muon_15_phi_vp[7]->Fill(phi);
                if(isAside){
                  th1_l1muon_15_phi_vp[18]->Fill(phi);
                  th1_l1muon_15_phi_vp[19]->Fill(phi);
                }else{
                  th1_l1muon_15_phi_vp[30]->Fill(phi);
                  th1_l1muon_15_phi_vp[31]->Fill(phi);
                }
              }
              if(isEndcap){
                th1_l1muon_15_pt_vp[10]->Fill(muon_pt);
                th1_l1muon_15_pt_vp[11]->Fill(muon_pt);
                th1_l1muon_15_phi_vp[10]->Fill(phi);
                th1_l1muon_15_phi_vp[11]->Fill(phi);
                if(isAside){
                  th1_l1muon_15_phi_vp[22]->Fill(phi);
                  th1_l1muon_15_phi_vp[23]->Fill(phi);
                }
                else{
                  th1_l1muon_15_phi_vp[34]->Fill(phi);
                  th1_l1muon_15_phi_vp[35]->Fill(phi);
                }
              }
              break;
            }

            case 2 /* loose */ :{
              th1_l1muon_15_eta_vp[3]->Fill(eta);
              th1_l1muon_15_phi_vp[3]->Fill(phi);
              if(isAside){  th1_l1muon_15_phi_vp[15]->Fill(phi);}
              else{         th1_l1muon_15_phi_vp[27]->Fill(phi);}
              th1_l1muon_15_pt_vp[3]->Fill(muon_pt);
              if(isBarrel){
                th1_l1muon_15_pt_vp[7]->Fill(muon_pt);
                th1_l1muon_15_phi_vp[7]->Fill(phi);
                if(isAside){  th1_l1muon_15_phi_vp[19]->Fill(phi);}
                else{         th1_l1muon_15_phi_vp[31]->Fill(phi);}
              }
              if(isEndcap){
                th1_l1muon_15_pt_vp[11]->Fill(muon_pt);
                th1_l1muon_15_phi_vp[11]->Fill(phi);
                if(isAside){  th1_l1muon_15_phi_vp[23]->Fill(phi);}
                else{         th1_l1muon_15_phi_vp[35]->Fill(phi);}
              }
              break;
            }
          } // end of quality switch
        } // end of L1MU15

        // for L1MU20
        if (pt_thr >= 6){
          th1_l1muon_20_eta_vp[0]->Fill(eta);
          th1_l1muon_20_phi_vp[0]->Fill(phi);
          th1_l1muon_20_pt_vp[0]->Fill(muon_pt);
          if(isBarrel){
            th1_l1muon_20_phi_vp[4]->Fill(phi);
            th1_l1muon_20_pt_vp[4]->Fill(muon_pt);
            if(isAside){  th1_l1muon_20_phi_vp[16]->Fill(phi);}
            else{         th1_l1muon_20_phi_vp[28]->Fill(phi);}
          }
          if(isEndcap){
            th1_l1muon_20_phi_vp[8]->Fill(phi);
            th1_l1muon_20_pt_vp[8]->Fill(muon_pt);
            if(isAside){  th1_l1muon_20_phi_vp[20]->Fill(phi);}
            else{         th1_l1muon_20_phi_vp[32]->Fill(phi);}
          }
          if(isAside){  th1_l1muon_20_phi_vp[12]->Fill(phi);}
          else{         th1_l1muon_20_phi_vp[24]->Fill(phi);}

          // Quality criteria
          switch (quality) {
            case 0 /* tight */ :{
              th1_l1muon_20_eta_vp[1]->Fill(eta);
              th1_l1muon_20_eta_vp[2]->Fill(eta);
              th1_l1muon_20_eta_vp[3]->Fill(eta);
              th1_l1muon_20_phi_vp[1]->Fill(phi);
              th1_l1muon_20_phi_vp[2]->Fill(phi);
              th1_l1muon_20_phi_vp[3]->Fill(phi);
              if(isAside){
                th1_l1muon_20_phi_vp[13]->Fill(phi);
                th1_l1muon_20_phi_vp[14]->Fill(phi);
                th1_l1muon_20_phi_vp[15]->Fill(phi);
              }
              else{
                th1_l1muon_20_phi_vp[25]->Fill(phi);
                th1_l1muon_20_phi_vp[26]->Fill(phi);
                th1_l1muon_20_phi_vp[27]->Fill(phi);
              }
              th1_l1muon_20_pt_vp[1]->Fill(muon_pt);
              th1_l1muon_20_pt_vp[2]->Fill(muon_pt);
              th1_l1muon_20_pt_vp[3]->Fill(muon_pt);
              if(isBarrel){
                th1_l1muon_20_pt_vp[5]->Fill(muon_pt);
                th1_l1muon_20_pt_vp[6]->Fill(muon_pt);
                th1_l1muon_20_pt_vp[7]->Fill(muon_pt);
                th1_l1muon_20_phi_vp[5]->Fill(phi);
                th1_l1muon_20_phi_vp[6]->Fill(phi);
                th1_l1muon_20_phi_vp[7]->Fill(phi);
                if(isAside){
                  th1_l1muon_20_phi_vp[17]->Fill(phi);
                  th1_l1muon_20_phi_vp[18]->Fill(phi);
                  th1_l1muon_20_phi_vp[19]->Fill(phi);
                }
                else{
                  th1_l1muon_20_phi_vp[29]->Fill(phi);
                  th1_l1muon_20_phi_vp[30]->Fill(phi);
                  th1_l1muon_20_phi_vp[31]->Fill(phi);
                }
              }
              if(isEndcap){
                th1_l1muon_20_pt_vp[9]->Fill(muon_pt);
                th1_l1muon_20_pt_vp[10]->Fill(muon_pt);
                th1_l1muon_20_pt_vp[11]->Fill(muon_pt);
                th1_l1muon_20_phi_vp[9]->Fill(phi);
                th1_l1muon_20_phi_vp[10]->Fill(phi);
                th1_l1muon_20_phi_vp[11]->Fill(phi);
                if(isAside){
                  th1_l1muon_20_phi_vp[21]->Fill(phi);
                  th1_l1muon_20_phi_vp[22]->Fill(phi);
                  th1_l1muon_20_phi_vp[23]->Fill(phi);
                }
                else{
                  th1_l1muon_20_phi_vp[33]->Fill(phi);
                  th1_l1muon_20_phi_vp[34]->Fill(phi);
                  th1_l1muon_20_phi_vp[35]->Fill(phi);
                }
              }
              break;
            }
            case 1 /* medium */:{
              th1_l1muon_20_eta_vp[2]->Fill(eta);
              th1_l1muon_20_eta_vp[3]->Fill(eta);
              th1_l1muon_20_phi_vp[2]->Fill(phi);
              th1_l1muon_20_phi_vp[3]->Fill(phi);
              if(isAside){
                th1_l1muon_20_phi_vp[14]->Fill(phi);
                th1_l1muon_20_phi_vp[15]->Fill(phi);
              }
              else{
                th1_l1muon_20_phi_vp[26]->Fill(phi);
                th1_l1muon_20_phi_vp[27]->Fill(phi);
              }
              th1_l1muon_20_pt_vp[2]->Fill(muon_pt);
              th1_l1muon_20_pt_vp[3]->Fill(muon_pt);
              if(isBarrel){
                th1_l1muon_20_pt_vp[6]->Fill(muon_pt);
                th1_l1muon_20_pt_vp[7]->Fill(muon_pt);
                th1_l1muon_20_phi_vp[6]->Fill(phi);
                th1_l1muon_20_phi_vp[7]->Fill(phi);
                if(isAside){
                  th1_l1muon_20_phi_vp[18]->Fill(phi);
                  th1_l1muon_20_phi_vp[19]->Fill(phi);
                }
                else{
                  th1_l1muon_20_phi_vp[30]->Fill(phi);
                  th1_l1muon_20_phi_vp[31]->Fill(phi);
                }
              }
              if(isEndcap){
                th1_l1muon_20_pt_vp[10]->Fill(muon_pt);
                th1_l1muon_20_pt_vp[11]->Fill(muon_pt);
                th1_l1muon_20_phi_vp[10]->Fill(phi);
                th1_l1muon_20_phi_vp[11]->Fill(phi);
                if(isAside){
                  th1_l1muon_20_phi_vp[22]->Fill(phi);
                  th1_l1muon_20_phi_vp[23]->Fill(phi);
                }
                else{
                  th1_l1muon_20_phi_vp[34]->Fill(phi);
                  th1_l1muon_20_phi_vp[35]->Fill(phi);
                }
              }
              break;
            }
            case 2 /* loose */ :{
              th1_l1muon_20_eta_vp[3]->Fill(eta);
              th1_l1muon_20_phi_vp[3]->Fill(phi);
              if(isAside){  th1_l1muon_20_phi_vp[15]->Fill(phi);}
              else{         th1_l1muon_20_phi_vp[27]->Fill(phi);}
              th1_l1muon_20_pt_vp[3]->Fill(muon_pt);
              if(isBarrel){
                th1_l1muon_20_pt_vp[7]->Fill(muon_pt);
                th1_l1muon_20_phi_vp[7]->Fill(phi);
                if(isAside){  th1_l1muon_20_phi_vp[19]->Fill(phi);}
                else{         th1_l1muon_20_phi_vp[31]->Fill(phi);}
              }
              if(isEndcap){
                th1_l1muon_20_pt_vp[11]->Fill(muon_pt);
                th1_l1muon_20_phi_vp[11]->Fill(phi);
                if(isAside){  th1_l1muon_20_phi_vp[23]->Fill(phi);}
                else{         th1_l1muon_20_phi_vp[35]->Fill(phi);}
              }
              break;
            }
          } // end of quality switch
        } // end of L1MU20

      } // end of true positive

    } // end of muon loop
  } // end of L1 Trigger nloop

} // End of l1_trigger_loop function

void TMDB::loop_tmdb( bool ignoreAside ){

  Int_t mf_out_array[total_tile_sides][total_tile_modules][total_tmdb_mod_ch];
  Int_t mf_decision_cell_array[total_tile_sides][total_tile_modules];
  Int_t mf_decision_mod_array[total_tile_sides][total_tile_modules];
  for(Int_t side=0; side<total_tile_sides; side++) { // side loop
    for(Int_t module=0; module<total_tile_modules; module++) { // module loop
      mf_decision_cell_array[side][module]=0;
      mf_decision_mod_array[side][module]=0;
    }
  }

  for(Int_t tmdb_raw_n=0; tmdb_raw_n<TILE_murcv_raw_n; tmdb_raw_n++){ // raw loop
    Int_t module = TILE_murcv_raw_drawer->at(tmdb_raw_n);
    Int_t side = TILE_murcv_raw_ros->at(tmdb_raw_n); // 3 = Aside, 4 = Cside
    Int_t channel = TILE_murcv_raw_channel->at(tmdb_raw_n);
    Int_t mf_out = TILE_murcv_raw_count->at(tmdb_raw_n);
    Float_t energy = TILE_murcv_raw_energy->at(tmdb_raw_n);

    if(side > 4 || side < 3){  continue;}
    side = side - 3; // 0 = Aside, 1 = Cside

    mf_out_array[side][module][channel] = mf_out;
    th1_mf_out[side*total_tile_modules*total_tmdb_mod_ch+module*total_tmdb_mod_ch+channel]->Fill(mf_out);
    th1_mf_energy[side*total_tile_modules*total_tmdb_mod_ch+module*total_tmdb_mod_ch+channel]->Fill(energy);

    for(Int_t tmdb_trig_n=0; tmdb_trig_n<TILE_murcv_trig_n; tmdb_trig_n++){ // trigger loop
      Int_t module_trig = TILE_murcv_trig_mod->at(tmdb_trig_n);
      Int_t side_trig = TILE_murcv_trig_part->at(tmdb_trig_n); // 3 = Aside, 4 = Cside

      if(side_trig > 4 || side_trig < 3){  continue;}
      side_trig = side_trig - 3; // 0 = Aside, 1 = Cside
      if(side!=side_trig){    continue;}
      if(module!=module_trig){continue;}
      // cout << "module: " << module << " module_trig: " << module_trig << " side: " << side << " side_trig: " << side_trig << endl;
      if(TILE_murcv_trig_bit1->at(tmdb_trig_n)==1){
        mf_decision_mod_array[side_trig][module_trig] = 1;
      }
      if(TILE_murcv_trig_bit3->at(tmdb_trig_n)==1){
        mf_decision_cell_array[side_trig][module_trig] = 1;
      }
    } // end of TMDB Trigger loop
  } // end of TMDB Raw loop

  // Filling Module and Cell histograms
  for(Int_t side=0; side<total_tile_sides; side++) { // side loop
    for(Int_t module=0; module<total_tile_modules; module++) { // module loop
      Int_t mf_mod_sum = mf_out_array[side][module][0]+mf_out_array[side][module][1]+mf_out_array[side][module][2]+mf_out_array[side][module][3];
      Int_t mf_cell_sum = mf_out_array[side][module][2]+mf_out_array[side][module][3];
      th1_mf_out_mod[side*total_tile_modules+module]->Fill(mf_mod_sum);
      th1_mf_out_cell[side*total_tile_modules+module]->Fill(mf_cell_sum);

      if(mf_decision_mod_array[side][module] == 1){
        th1_mf_decision_mod[side*total_tile_modules+module]->Fill(mf_mod_sum);
      }
      if(mf_decision_cell_array[side][module] == 1){
        th1_mf_decision_cell[side*total_tile_modules+module]->Fill(mf_cell_sum);
      }

    } // end of module loop
  } // end of side loop

} // end of tmdb_loop function

void TMDB::loop_extra( bool ignoreAside ){

  // L1 Trigger loop
  for(Int_t nloop=0; nloop<trig_L1_mu_n; nloop++) {
    Float_t eta   = trig_L1_mu_eta->at(nloop);
    Float_t phi   = trig_L1_mu_phi->at(nloop);
    Int_t pt_thr  = trig_L1_mu_thrNumber->at(nloop);
    Int_t roi     = trig_L1_mu_RoINumber->at(nloop);
    Int_t ssc     = (roi+4)/8;
    Int_t sector  = (trig_L1_mu_sectorAddress->at(nloop)>>1) & (0x3F);
    Int_t tile_module = tilemodule_to_check(sector);
    // Int_t vetoed  = trig_L1_mu_vetoed->at(nloop);
    bool isAside  = (trig_L1_mu_hemisphere->at(nloop)==1) ? true: false;
    // bool isAside  = (eta>0.0) ? true: false;
    // bool isBarrel = (trig_L1_mu_source->at(nloop)==0) ? true: false;
    bool isEndcap = (trig_L1_mu_source->at(nloop)==1) ? true: false;
    // if(!(isEndcap && pt_thr >= 5 && ssc < 6)) { continue; }
    // if( fabs(eta) < 0.9 || 1.3 < fabs(eta) )  { continue; }
    // if(vetoed) {continue;}
    if( ignoreAside && isAside) {continue;} // test
    bool threshold = false;

    // Get data from TMDB
    for(Int_t tmdb_n=0; tmdb_n<TILE_murcv_trig_n; tmdb_n++){
      Int_t module = TILE_murcv_trig_mod->at(tmdb_n);
      if (tile_module < module || tile_module > module+1) {continue;} // tile module different
      if(TILE_murcv_trig_bit1->at(tmdb_n) || TILE_murcv_trig_bit3->at(tmdb_n)) {
        threshold = true;
        th1_tmdb->Fill(eta);
        th1_tmdb_phi[0]->Fill(phi);

        if(isAside){th1_tmdb_phi[1]->Fill(phi);}
        else{       th1_tmdb_phi[2]->Fill(phi);}
      }
    } // end of tmdb coincidence loop

    // offline muon associated
    for(Int_t muon=0; muon<mu_n; muon++) {
      Float_t muon_eta = mu_eta->at(muon);
      Float_t muon_phi = mu_phi->at(muon);
      Float_t muon_pt  = mu_pt->at(muon);
      Float_t dR = calc_dR(eta-muon_eta , phi-muon_phi);
      // bool isCombined = (mu_muonType->at(muon)==0)? true: false;
      // bool isCaloTag = (mu_muonType->at(muon)==3)? true: false;
      Int_t quality	  = mu_quality->at(muon);
      // if(!isCombined)   {continue;} // combined muon cut
      // if(dR> 0.1)       {continue;} // dR cut

      // true positive: it has an offline muon associated L1MU15
      if(dR < 0.1 && threshold){
        th1_tmdb_15_eta_vp[0]->Fill(muon_eta);
        th1_tmdb_15_phi_vp[0]->Fill(muon_phi);
        th1_tmdb_15_pt_vp[0]->Fill(muon_pt);
        if(isAside){th1_tmdb_15_phi_vp[4]->Fill(muon_phi);}
        else{       th1_tmdb_15_phi_vp[8]->Fill(muon_phi);}
        // Quality criteria
        switch (quality) {
          case 0 /* tight */ :{
            th1_tmdb_15_eta_vp[1]->Fill(muon_eta);
            th1_tmdb_15_eta_vp[2]->Fill(muon_eta);
            th1_tmdb_15_eta_vp[3]->Fill(muon_eta);
            th1_tmdb_15_phi_vp[1]->Fill(muon_phi);
            th1_tmdb_15_phi_vp[2]->Fill(muon_phi);
            th1_tmdb_15_phi_vp[3]->Fill(muon_phi);
            th1_tmdb_15_pt_vp[1]->Fill(muon_pt);
            th1_tmdb_15_pt_vp[2]->Fill(muon_pt);
            th1_tmdb_15_pt_vp[3]->Fill(muon_pt);
            if(isAside){
              th1_tmdb_15_phi_vp[5]->Fill(muon_phi);
              th1_tmdb_15_phi_vp[6]->Fill(muon_phi);
              th1_tmdb_15_phi_vp[7]->Fill(muon_phi);
            }else{
              th1_tmdb_15_phi_vp[9]->Fill(muon_phi);
              th1_tmdb_15_phi_vp[10]->Fill(muon_phi);
              th1_tmdb_15_phi_vp[11]->Fill(muon_phi);
            }
            break;
          }
          case 1 /* medium */:{
            th1_tmdb_15_eta_vp[2]->Fill(muon_eta);
            th1_tmdb_15_eta_vp[3]->Fill(muon_eta);
            th1_tmdb_15_phi_vp[2]->Fill(muon_phi);
            th1_tmdb_15_phi_vp[3]->Fill(muon_phi);
            th1_tmdb_15_pt_vp[2]->Fill(muon_pt);
            th1_tmdb_15_pt_vp[3]->Fill(muon_pt);
            if(isAside){
              th1_tmdb_15_phi_vp[6]->Fill(muon_phi);
              th1_tmdb_15_phi_vp[7]->Fill(muon_phi);
            }else{
              th1_tmdb_15_phi_vp[10]->Fill(muon_phi);
              th1_tmdb_15_phi_vp[11]->Fill(muon_phi);
            }
            break;
          }
          case 2 /* loose */ :{
            th1_tmdb_15_eta_vp[3]->Fill(muon_eta);
            th1_tmdb_15_phi_vp[3]->Fill(muon_phi);
            th1_tmdb_15_pt_vp[3]->Fill(muon_pt);
            if(isAside){th1_tmdb_15_phi_vp[7]->Fill(muon_phi);}
            else{       th1_tmdb_15_phi_vp[11]->Fill(muon_phi);}
            break;
          }
        } // end of quality switch
      } // end of true positive for L1MU15

      // true positive: it has an offline muon associated L1MU20
      if(dR < 0.1 && threshold && pt_thr>=6){
        th1_tmdb_20_eta_vp[0]->Fill(muon_eta);
        th1_tmdb_20_phi_vp[0]->Fill(muon_phi);
        th1_tmdb_20_pt_vp[0]->Fill(muon_pt);
        if(isAside){th1_tmdb_20_phi_vp[4]->Fill(muon_phi);}
        else{       th1_tmdb_20_phi_vp[8]->Fill(muon_phi);}
        // Quality criteria
        switch (quality) {
          case 0 /* tight */ :{
            th1_tmdb_20_eta_vp[1]->Fill(muon_eta);
            th1_tmdb_20_eta_vp[2]->Fill(muon_eta);
            th1_tmdb_20_eta_vp[3]->Fill(muon_eta);
            th1_tmdb_20_phi_vp[1]->Fill(muon_phi);
            th1_tmdb_20_phi_vp[2]->Fill(muon_phi);
            th1_tmdb_20_phi_vp[3]->Fill(muon_phi);
            th1_tmdb_20_pt_vp[1]->Fill(muon_pt);
            th1_tmdb_20_pt_vp[2]->Fill(muon_pt);
            th1_tmdb_20_pt_vp[3]->Fill(muon_pt);
            if(isAside){
              th1_tmdb_20_phi_vp[5]->Fill(muon_phi);
              th1_tmdb_20_phi_vp[6]->Fill(muon_phi);
              th1_tmdb_20_phi_vp[7]->Fill(muon_phi);
            }else{
              th1_tmdb_20_phi_vp[9]->Fill(muon_phi);
              th1_tmdb_20_phi_vp[10]->Fill(muon_phi);
              th1_tmdb_20_phi_vp[11]->Fill(muon_phi);
            }
            break;
          }
          case 1 /* medium */:{
            th1_tmdb_20_eta_vp[2]->Fill(muon_eta);
            th1_tmdb_20_eta_vp[3]->Fill(muon_eta);
            th1_tmdb_20_phi_vp[2]->Fill(muon_phi);
            th1_tmdb_20_phi_vp[3]->Fill(muon_phi);
            th1_tmdb_20_pt_vp[2]->Fill(muon_pt);
            th1_tmdb_20_pt_vp[3]->Fill(muon_pt);
            if(isAside){
              th1_tmdb_20_phi_vp[6]->Fill(muon_phi);
              th1_tmdb_20_phi_vp[7]->Fill(muon_phi);
            }else{
              th1_tmdb_20_phi_vp[10]->Fill(muon_phi);
              th1_tmdb_20_phi_vp[11]->Fill(muon_phi);
            }
            break;
          }
          case 2 /* loose */ :{
            th1_tmdb_20_eta_vp[3]->Fill(muon_eta);
            th1_tmdb_20_phi_vp[3]->Fill(muon_phi);
            th1_tmdb_20_pt_vp[3]->Fill(muon_pt);
            if(isAside){th1_tmdb_20_phi_vp[7]->Fill(muon_phi);}
            else{       th1_tmdb_20_phi_vp[11]->Fill(muon_phi);}
            break;
          }
        } // end of quality switch
      } // end of true positive for L1MU20

    } // end of offline loop

  } // end of L1Trigger loop
}

// ###############################################################################
// Int_t TMDB::loop_save_muons(string fout)
// Description: Save muons in text files
// Requirements: folder name
// ###############################################################################
void TMDB::loop_save_muons(string fout){

  for(Int_t tile=0; tile<TILE_cell_n; tile++){ // trigger loop
    // Float_t energy    = TILE_cell_E->at(tile); // Energy of cell
    Float_t ene1      = TILE_cell_ene1->at(tile); // energy of one cell PMT
    Float_t eta       = TILE_cell_eta->at(tile); // energy of one cell PMT
    Int_t partition   = TILE_cell_partition->at(tile); // EBA:3 or EBC:4
    Int_t section     = TILE_cell_section->at(tile);
    Int_t module_tile = TILE_cell_module->at(tile);
    Int_t cell        = TILE_cell_sample->at(tile);

    if(partition<3){continue;}  // cut for only EBA:3 and EBC:4
    if(section!=2){continue;}   // cut for only extended barrel
    if(cell!=2){continue;}      // cut for only Dlayer cells
    if(ene1 < 2000 || ene1 > 5000){continue;} // cut energy to be a muon
    bool isD5 = (abs(eta)<1.2) ? true: false;

    for(Int_t tmdb=0; tmdb<TILE_murcv_digit_n; tmdb++){ // raw loop
      Int_t side        = TILE_murcv_digit_ros->at(tmdb); // Aside:3 or Cside:4
      Int_t module_tmdb = TILE_murcv_digit_drawer->at(tmdb); // module number
      Int_t channel     = TILE_murcv_digit_channel->at(tmdb); // 0=D5L, 1=D5R, 2=D6L , 3=D6R
      if(module_tile!=module_tmdb){continue;} // cut for same module
      if(partition!=side){continue;}          // cut for same side
      if( isD5 && channel>=2){continue;}      // cut for same cell D5
      if(!isD5 && channel<=1){continue;}      // cut for same cell D6

      // TMDB Muon samples
      vector<float> sampleTMDB = TILE_murcv_digit_sampleVec->at(tmdb);
      string sidestr;if(side==3){sidestr = "sideA";}else{sidestr="sideC";}
      string outputfile = fout + "/muons/" + sidestr + "/mod_" + std::to_string(module_tmdb) +"_ch_" + std::to_string(channel);
      ofstream of; of.open(outputfile.c_str(), std::ofstream::app); // open file
      for (Int_t sample=0; sample < total_tmdb_samples; sample++){ // Loop to write sampleTMDB converting Uchar_t* to Double_t*
        of << sampleTMDB[sample];
        if(sample < 6){of << std::setw(5);}
      } // end of sample loop
      of << std::endl; // end of line
      of.close(); // close stream
    } // end of digit loop
  } // end of Tile Cell loop

} // End of loop_save_muons

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
  // for (Int_t i=0; i<5; i++){th1_offline[i]->Write();}
  for (Int_t i=0; i<9; i++){th1_offline_eta[i]->Write();}
  for (Int_t i=0; i<12; i++){th1_offline_phi[i]->Write();}
  for (Int_t i=0; i<16; i++){th1_offline_pt[i]->Write();}

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
  // for (Int_t i=0; i<4; i++){  th1_l1muon_15_eta_vn[i]->Write();}
  // for (Int_t i=0; i<4; i++){  th1_l1muon_15_eta_fp[i]->Write();}
  // for (Int_t i=0; i<4; i++){  th1_l1muon_15_eta_fn[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_15_phi[i]->Write();}
  for (Int_t i=0; i<36; i++){ th1_l1muon_15_phi_vp[i]->Write();}
  // for (Int_t i=0; i<36; i++){ th1_l1muon_15_phi_vn[i]->Write();}
  // for (Int_t i=0; i<36; i++){ th1_l1muon_15_phi_fp[i]->Write();}
  // for (Int_t i=0; i<36; i++){ th1_l1muon_15_phi_fn[i]->Write();}
  for (Int_t i=0; i<3; i++){  th1_l1muon_15_pt[i]->Write();}
  for (Int_t i=0; i<12; i++){ th1_l1muon_15_pt_vp[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_l1muon_15_pt_vn[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_l1muon_15_pt_fp[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_l1muon_15_pt_fn[i]->Write();}
  th1_l1muon_20_eta->Write();
  for (Int_t i=0; i<4; i++){  th1_l1muon_20_eta_vp[i]->Write();}
  // for (Int_t i=0; i<4; i++){  th1_l1muon_20_eta_vn[i]->Write();}
  // for (Int_t i=0; i<4; i++){  th1_l1muon_20_eta_fp[i]->Write();}
  // for (Int_t i=0; i<4; i++){  th1_l1muon_20_eta_fn[i]->Write();}
  for (Int_t i=0; i<9; i++){  th1_l1muon_20_phi[i]->Write();}
  for (Int_t i=0; i<36; i++){ th1_l1muon_20_phi_vp[i]->Write();}
  // for (Int_t i=0; i<36; i++){ th1_l1muon_20_phi_vn[i]->Write();}
  // for (Int_t i=0; i<36; i++){ th1_l1muon_20_phi_fp[i]->Write();}
  // for (Int_t i=0; i<36; i++){ th1_l1muon_20_phi_fn[i]->Write();}
  for (Int_t i=0; i<3; i++){  th1_l1muon_20_pt[i]->Write();}
  for (Int_t i=0; i<12; i++){ th1_l1muon_20_pt_vp[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_l1muon_20_pt_vn[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_l1muon_20_pt_fp[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_l1muon_20_pt_fn[i]->Write();}

  // Writing TMDB histograms
  // th1_tmdb->Write();
  th1_tmdb->Write();
  th1_tmdb_eta->Write();
  for (Int_t i=0; i<3; i++){  th1_tmdb_phi[i]->Write();}
  for (Int_t i=0; i<5; i++){  th1_tmdb_15_eta_vp[i]->Write();}
  // for (Int_t i=0; i<5; i++){  th1_tmdb_15_eta_vn[i]->Write();}
  // for (Int_t i=0; i<5; i++){  th1_tmdb_15_eta_fp[i]->Write();}
  // for (Int_t i=0; i<5; i++){  th1_tmdb_15_eta_fn[i]->Write();}
  // for (Int_t i=0; i<3; i++){  th1_tmdb_15_phi[i]->Write();}
  for (Int_t i=0; i<12; i++){ th1_tmdb_15_phi_vp[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_tmdb_15_phi_vn[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_tmdb_15_phi_fp[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_tmdb_15_phi_fn[i]->Write();}
  for (Int_t i=0; i<5; i++){ th1_tmdb_15_pt_vp[i]->Write();}

  // th1_tmdb_20->Write();
  // th1_tmdb_20_eta->Write();
  for (Int_t i=0; i<5; i++){  th1_tmdb_20_eta_vp[i]->Write();}
  // for (Int_t i=0; i<5; i++){  th1_tmdb_20_eta_vn[i]->Write();}
  // for (Int_t i=0; i<5; i++){  th1_tmdb_20_eta_fp[i]->Write();}
  // for (Int_t i=0; i<5; i++){  th1_tmdb_20_eta_fn[i]->Write();}
  // for (Int_t i=0; i<3; i++){  th1_tmdb_20_phi[i]->Write();}
  for (Int_t i=0; i<12; i++){ th1_tmdb_20_phi_vp[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_tmdb_20_phi_vn[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_tmdb_20_phi_fp[i]->Write();}
  // for (Int_t i=0; i<12; i++){ th1_tmdb_20_phi_fn[i]->Write();}
  for (Int_t i=0; i<5; i++){ th1_tmdb_20_pt_vp[i]->Write();}

  // MF output
  for (Int_t i=0; i<512; i++){ th1_mf_out[i]->Write();}
  for (Int_t i=0; i<512; i++){ th1_mf_energy[i]->Write();}
  for (Int_t i=0; i<128; i++){ th1_mf_out_cell[i]->Write();}
  for (Int_t i=0; i<128; i++){ th1_mf_out_mod[i]->Write();}
  for (Int_t i=0; i<128; i++){ th1_mf_decision_cell[i]->Write();}
  for (Int_t i=0; i<128; i++){ th1_mf_decision_mod[i]->Write();}


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
