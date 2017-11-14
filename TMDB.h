//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 30 14:55:32 2015 by ROOT version 6.02/12
// from TTree TMDB/TMDB
// found on file: L1TGCNtuple.root
//////////////////////////////////////////////////////////

#ifndef TMDB_h
#define TMDB_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TPrincipal.h>
#include <TRobustEstimator.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDEigen.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TLegend.h>

// Header file for the classes stored in the TTree if any.
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <bitset>
#include <vector>
#include <string>
#include <algorithm>
#include <vector>


#define total_tile_modules		64 //64
#define total_tile_sides 		2
#define total_tmdb_mod_ch		4 //4
#define total_tmdb_samples		7
#define total_tmdb_coeff		7
#define cov_matrix_size			7
#define total_fit_const			2

using namespace std;

class TMDB {
private :


public :
   TTree          *fChainPhysics;   //!pointer to the analyzed Physics TTree or TChain
   Int_t           fCurrentPhysics; //!current Physics Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types for h2000
   UChar_t         sampleTMDB[4][64][8][7];	// sampleTMDB[side][module][channel][sample] // side: 0->LBA, 1->LBC, 2->EBA, 3->EBC
//	Float_t         eDsp[4][64][48];				// for h2000 physics runs

   // Declaration of leaf types for L1TGCNtuple
   Int_t           RunNumber;
   Int_t           EventNumber;
   Int_t           timeStamp;
   Int_t           timeStampNSOffset;
   Int_t           lbn;
   Int_t           bcid;
   Int_t           detmask0;
   Int_t           detmask1;
   Float_t         actualIntPerXing;
   Float_t         averageIntPerXing;
   Int_t           pixelFlags;
   Int_t           sctFlags;
   Int_t           trtFlags;
   Int_t           larFlags;
   Int_t           tileFlags;
   Int_t           muonFlags;
   Int_t           fwdFlags;
   Int_t           coreFlags;
   Int_t           pixelError;
   Int_t           sctError;
   Int_t           trtError;
   Int_t           larError;
   Int_t           tileError;
   Int_t           muonError;
   Int_t           fwdError;
   Int_t           coreError;
   Int_t           trig_L1_mu_n;
   vector<float>   *trig_L1_mu_eta;
   vector<float>   *trig_L1_mu_phi;
   vector<string>  *trig_L1_mu_thrName;
   vector<short>   *trig_L1_mu_thrNumber;
   vector<short>   *trig_L1_mu_RoINumber;
   vector<short>   *trig_L1_mu_sectorAddress;
   vector<int>     *trig_L1_mu_firstCandidate;
   vector<int>     *trig_L1_mu_moreCandInRoI;
   vector<int>     *trig_L1_mu_moreCandInSector;
   vector<short>   *trig_L1_mu_source;
   vector<short>   *trig_L1_mu_hemisphere;
   vector<short>   *trig_L1_mu_charge;
   vector<int>     *trig_L1_mu_vetoed;
   Int_t           mu_n;
   vector<float>   *mu_pt;
   vector<float>   *mu_eta;
   vector<float>   *mu_phi;
   vector<float>   *mu_m;
   vector<int>     *mu_charge;
   vector<int>     *mu_author;
   vector<unsigned short> *mu_allAuthors;
   vector<int>     *mu_muonType;
   vector<float>   *mu_etcone20;
   vector<float>   *mu_etcone30;
   vector<float>   *mu_etcone40;
   vector<float>   *mu_ptcone20;
   vector<float>   *mu_ptcone30;
   vector<float>   *mu_ptcone40;
   vector<int>     *mu_isPassedMCP;
   vector<int>     *mu_quality;
   vector<float>   *mu_trackfitchi2;
   vector<float>   *mu_trackfitndof;
   vector<float>   *mu_msInnerMatchChi2;
   vector<float>   *mu_msOuterMatchChi2;
   vector<int>     *mu_msInnerMatchDOF;
   vector<int>     *mu_msOuterMatchDOF;
   vector<int>     *mu_nOutliersOnTrack;
   vector<int>     *mu_nBLHits;
   vector<int>     *mu_nPixHits;
   vector<int>     *mu_nSCTHits;
   vector<int>     *mu_nTRTHits;
   vector<int>     *mu_nTRTHighTHits;
   vector<int>     *mu_nBLSharedHits;
   vector<int>     *mu_nPixSharedHits;
   vector<int>     *mu_nPixHoles;
   vector<int>     *mu_nSCTSharedHits;
   vector<int>     *mu_nSCTHoles;
   vector<int>     *mu_nTRTOutliers;
   vector<int>     *mu_nTRTHighTOutliers;
   vector<int>     *mu_nGangedPixels;
   vector<int>     *mu_nPixelDeadSensors;
   vector<int>     *mu_nSCTDeadSensors;
   vector<int>     *mu_nTRTDeadStraws;
   vector<int>     *mu_expectBLayerHit;
   vector<int>     *mu_nPrecisionLayers;
   vector<int>     *mu_nPrecisionHoleLayers;
   vector<int>     *mu_nPhiLayers;
   vector<int>     *mu_nPhiHoleLayers;
   vector<int>     *mu_nTrigEtaLayers;
   vector<int>     *mu_nTrigEtaHoleLayers;
   vector<int>     *mu_primarySector;
   vector<int>     *mu_secondarySector;
   vector<int>     *mu_nInnerSmallHits;
   vector<int>     *mu_nInnerLargeHits;
   vector<int>     *mu_nMiddleSmallHits;
   vector<int>     *mu_nMiddleLargeHits;
   vector<int>     *mu_nOuterSmallHits;
   vector<int>     *mu_nOuterLargeHits;
   vector<int>     *mu_nExtendedSmallHits;
   vector<int>     *mu_nExtendedLargeHits;
   vector<int>     *mu_nInnerSmallHoles;
   vector<int>     *mu_nInnerLargeHoles;
   vector<int>     *mu_nMiddleSmallHoles;
   vector<int>     *mu_nMiddleLargeHoles;
   vector<int>     *mu_nOuterSmallHoles;
   vector<int>     *mu_nOuterLargeHoles;
   vector<int>     *mu_nExtendedSmallHoles;
   vector<int>     *mu_nExtendedLargeHoles;
   vector<int>     *mu_nPhiLayer1Hits;
   vector<int>     *mu_nPhiLayer2Hits;
   vector<int>     *mu_nPhiLayer3Hits;
   vector<int>     *mu_nPhiLayer4Hits;
   vector<int>     *mu_nEtaLayer1Hits;
   vector<int>     *mu_nEtaLayer2Hits;
   vector<int>     *mu_nEtaLayer3Hits;
   vector<int>     *mu_nEtaLayer4Hits;
   vector<int>     *mu_nPhiLayer1Holes;
   vector<int>     *mu_nPhiLayer2Holes;
   vector<int>     *mu_nPhiLayer3Holes;
   vector<int>     *mu_nPhiLayer4Holes;
   vector<int>     *mu_nEtaLayer1Holes;
   vector<int>     *mu_nEtaLayer2Holes;
   vector<int>     *mu_nEtaLayer3Holes;
   vector<int>     *mu_nEtaLayer4Holes;
   vector<float>   *mu_cb_d0;
   vector<float>   *mu_cb_z0;
   vector<float>   *mu_cb_phi0;
   vector<float>   *mu_cb_theta;
   vector<float>   *mu_cb_qOverP;
   vector<float>   *mu_cb_vx;
   vector<float>   *mu_cb_vy;
   vector<float>   *mu_cb_vz;
   Int_t           museg_n;
   vector<float>   *museg_x;
   vector<float>   *museg_y;
   vector<float>   *museg_z;
   vector<float>   *museg_px;
   vector<float>   *museg_py;
   vector<float>   *museg_pz;
   vector<float>   *museg_t0;
   vector<float>   *museg_t0error;
   vector<float>   *museg_chi2;
   vector<float>   *museg_ndof;
   vector<int>     *museg_sector;
   vector<int>     *museg_stationName;
   vector<int>     *museg_stationEta;
   vector<int>     *museg_author;
   Int_t           ext_mu_bias_n;
   vector<int>     *ext_mu_bias_type;
   vector<int>     *ext_mu_bias_index;
   vector<int>     *ext_mu_bias_size;
   vector<vector<int> > *ext_mu_bias_targetVec;
   vector<vector<float> > *ext_mu_bias_targetDistanceVec;
   vector<vector<float> > *ext_mu_bias_targetEtaVec;
   vector<vector<float> > *ext_mu_bias_targetPhiVec;
   vector<vector<float> > *ext_mu_bias_targetDeltaEtaVec;
   vector<vector<float> > *ext_mu_bias_targetDeltaPhiVec;
   vector<vector<float> > *ext_mu_bias_targetPxVec;
   vector<vector<float> > *ext_mu_bias_targetPyVec;
   vector<vector<float> > *ext_mu_bias_targetPzVec;
   Int_t           ext_mu_ubias_n;
   vector<int>     *ext_mu_ubias_type;
   vector<int>     *ext_mu_ubias_index;
   vector<int>     *ext_mu_ubias_size;
   vector<vector<int> > *ext_mu_ubias_targetVec;
   vector<vector<float> > *ext_mu_ubias_targetDistanceVec;
   vector<vector<float> > *ext_mu_ubias_targetEtaVec;
   vector<vector<float> > *ext_mu_ubias_targetPhiVec;
   vector<vector<float> > *ext_mu_ubias_targetDeltaEtaVec;
   vector<vector<float> > *ext_mu_ubias_targetDeltaPhiVec;
   vector<vector<float> > *ext_mu_ubias_targetPxVec;
   vector<vector<float> > *ext_mu_ubias_targetPyVec;
   vector<vector<float> > *ext_mu_ubias_targetPzVec;
   Int_t           trigger_info_n;
   vector<string>  *trigger_info_chain;
   vector<int>     *trigger_info_isPassed;
   vector<int>     *trigger_info_nTracks;
   vector<vector<int> > *trigger_info_typeVec;
   vector<vector<float> > *trigger_info_ptVec;
   vector<vector<float> > *trigger_info_etaVec;
   vector<vector<float> > *trigger_info_phiVec;
   Int_t           vxp_n;
   vector<float>   *vxp_x;
   vector<float>   *vxp_y;
   vector<float>   *vxp_z;
   vector<float>   *vxp_cov_x;
   vector<float>   *vxp_cov_y;
   vector<float>   *vxp_cov_z;
   vector<float>   *vxp_cov_xy;
   vector<float>   *vxp_cov_xz;
   vector<float>   *vxp_cov_yz;
   vector<float>   *vxp_chi2;
   vector<int>     *vxp_ndof;
   vector<int>     *vxp_nTracks;
   vector<int>     *vxp_type;
   Int_t           TGC_prd_n;
   vector<float>   *TGC_prd_x;
   vector<float>   *TGC_prd_y;
   vector<float>   *TGC_prd_z;
   vector<float>   *TGC_prd_shortWidth;
   vector<float>   *TGC_prd_longWidth;
   vector<float>   *TGC_prd_length;
   vector<int>     *TGC_prd_isStrip;
   vector<int>     *TGC_prd_gasGap;
   vector<int>     *TGC_prd_channel;
   vector<int>     *TGC_prd_eta;
   vector<int>     *TGC_prd_phi;
   vector<int>     *TGC_prd_station;
   vector<int>     *TGC_prd_bunch;
   Int_t           MDT_prd_n;
   vector<float>   *MDT_prd_x;
   vector<float>   *MDT_prd_y;
   vector<float>   *MDT_prd_z;
   vector<int>     *MDT_prd_adc;
   vector<int>     *MDT_prd_tdc;
   vector<int>     *MDT_prd_status;
   vector<float>   *MDT_prd_drift_radius;
   vector<float>   *MDT_prd_drift_radius_error;
   Int_t           TGC_coin_n;
   vector<float>   *TGC_coin_x_In;
   vector<float>   *TGC_coin_y_In;
   vector<float>   *TGC_coin_z_In;
   vector<float>   *TGC_coin_x_Out;
   vector<float>   *TGC_coin_y_Out;
   vector<float>   *TGC_coin_z_Out;
   vector<float>   *TGC_coin_width_In;
   vector<float>   *TGC_coin_width_Out;
   vector<float>   *TGC_coin_width_R;
   vector<float>   *TGC_coin_width_Phi;
   vector<int>     *TGC_coin_isAside;
   vector<int>     *TGC_coin_isForward;
   vector<int>     *TGC_coin_isStrip;
   vector<int>     *TGC_coin_isInner;
   vector<int>     *TGC_coin_isPositiveDeltaR;
   vector<int>     *TGC_coin_type;
   vector<int>     *TGC_coin_trackletId;
   vector<int>     *TGC_coin_trackletIdStrip;
   vector<int>     *TGC_coin_phi;
   vector<int>     *TGC_coin_roi;
   vector<int>     *TGC_coin_pt;
   vector<int>     *TGC_coin_delta;
   vector<int>     *TGC_coin_sub;
   vector<int>     *TGC_coin_veto;
   vector<int>     *TGC_coin_bunch;
   vector<int>     *TGC_coin_inner;
   Int_t           TILE_murcv_trig_n;
   vector<int>     *TILE_murcv_trig_mod;
   vector<int>     *TILE_murcv_trig_part;
   vector<bool>    *TILE_murcv_trig_bit0;
   vector<bool>    *TILE_murcv_trig_bit1;
   vector<bool>    *TILE_murcv_trig_bit2;
   vector<bool>    *TILE_murcv_trig_bit3;
   Int_t           TILE_murcv_raw_n;
   vector<float>   *TILE_murcv_raw_count;
   vector<float>   *TILE_murcv_raw_energy;
   vector<int>     *TILE_murcv_raw_ros;
   vector<int>     *TILE_murcv_raw_drawer;
   vector<int>     *TILE_murcv_raw_channel;
   Int_t           TILE_murcv_digit_n;
   vector<int>     *TILE_murcv_digit_nSamples;
   vector<int>     *TILE_murcv_digit_ros;
   vector<int>     *TILE_murcv_digit_drawer;
   vector<int>     *TILE_murcv_digit_channel;
   vector<vector<float> > *TILE_murcv_digit_sampleVec;
   Int_t           TGC_hierarchy_n;
   vector<int>     *TGC_hierarchy_index;
   vector<int>     *TGC_hierarchy_dR_hiPt;
   vector<int>     *TGC_hierarchy_dPhi_hiPt;
   vector<int>     *TGC_hierarchy_dR_tracklet;
   vector<int>     *TGC_hierarchy_dPhi_tracklet;
   vector<int>     *TGC_hierarchy_isChamberBoundary;
   vector<unsigned int> *muctpi_candidateMultiplicities;
   Int_t           muctpi_nDataWords;
   vector<unsigned int> *muctpi_dataWords;
   vector<float>   *muctpi_dw_eta;
   vector<float>   *muctpi_dw_phi;
   vector<short>   *muctpi_dw_source;
   vector<short>   *muctpi_dw_hemisphere;
   vector<short>   *muctpi_dw_bcid;
   vector<short>   *muctpi_dw_sectorID;
   vector<short>   *muctpi_dw_thrNumber;
   vector<short>   *muctpi_dw_roi;
   vector<short>   *muctpi_dw_veto;
   vector<short>   *muctpi_dw_firstCandidate;
   vector<short>   *muctpi_dw_moreCandInRoI;
   vector<short>   *muctpi_dw_moreCandInSector;
   vector<short>   *muctpi_dw_charge;
   vector<short>   *muctpi_dw_candidateVetoed;
  //  Int_t           TILE_cell_n;
  //  vector<float>   *TILE_cell_E;
  //  vector<float>   *TILE_cell_Et;
  //  vector<float>   *TILE_cell_eta;
  //  vector<float>   *TILE_cell_phi;
  //  vector<float>   *TILE_cell_sinTh;
  //  vector<float>   *TILE_cell_cosTh;
  //  vector<float>   *TILE_cell_cotTh;
  //  vector<float>   *TILE_cell_x;
  //  vector<float>   *TILE_cell_y;
  //  vector<float>   *TILE_cell_z;
  //  vector<int>     *TILE_cell_badcell;
  //  vector<int>     *TILE_cell_partition;
  //  vector<int>     *TILE_cell_section;
  //  vector<int>     *TILE_cell_side;
  //  vector<int>     *TILE_cell_module;
  //  vector<int>     *TILE_cell_tower;
  //  vector<int>     *TILE_cell_sample;

   // List of branches for L1TGCNtuple
   TBranch        *b_runNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_timeStamp;   //!
   TBranch        *b_timeStampNSOffset;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_detmask0;   //!
   TBranch        *b_detmask1;   //!
   TBranch        *b_actualIntPerXing;   //!
   TBranch        *b_averageIntPerXing;   //!
   TBranch        *b_pixelFlags;   //!
   TBranch        *b_sctFlags;   //!
   TBranch        *b_trtFlags;   //!
   TBranch        *b_larFlags;   //!
   TBranch        *b_tileFlags;   //!
   TBranch        *b_muonFlags;   //!
   TBranch        *b_fwdFlags;   //!
   TBranch        *b_coreFlags;   //!
   TBranch        *b_pixelError;   //!
   TBranch        *b_sctError;   //!
   TBranch        *b_trtError;   //!
   TBranch        *b_larError;   //!
   TBranch        *b_tileError;   //!
   TBranch        *b_muonError;   //!
   TBranch        *b_fwdError;   //!
   TBranch        *b_coreError;   //!
   TBranch        *b_trig_L1_mu_n;   //!
   TBranch        *b_trig_L1_mu_eta;   //!
   TBranch        *b_trig_L1_mu_phi;   //!
   TBranch        *b_trig_L1_mu_thrName;   //!
   TBranch        *b_trig_L1_mu_thrNumber;   //!
   TBranch        *b_trig_L1_mu_RoINumber;   //!
   TBranch        *b_trig_L1_mu_sectorAddress;   //!
   TBranch        *b_trig_L1_mu_firstCandidate;   //!
   TBranch        *b_trig_L1_mu_moreCandInRoI;   //!
   TBranch        *b_trig_L1_mu_moreCandInSector;   //!
   TBranch        *b_trig_L1_mu_source;   //!
   TBranch        *b_trig_L1_mu_hemisphere;   //!
   TBranch        *b_trig_L1_mu_charge;   //!
   TBranch        *b_trig_L1_mu_vetoed;   //!
   TBranch        *b_mu_n;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_m;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_author;   //!
   TBranch        *b_mu_allAuthors;   //!
   TBranch        *b_mu_muonType;   //!
   TBranch        *b_mu_etcone20;   //!
   TBranch        *b_mu_etcone30;   //!
   TBranch        *b_mu_etcone40;   //!
   TBranch        *b_mu_ptcone20;   //!
   TBranch        *b_mu_ptcone30;   //!
   TBranch        *b_mu_ptcone40;   //!
   TBranch        *b_mu_trackfitchi2;   //!
   TBranch        *b_mu_trackfitndof;   //!
   TBranch        *b_mu_isPassedMCP;   //!
   TBranch        *b_mu_quality;   //!
   TBranch        *b_mu_msInnerMatchChi2;   //!
   TBranch        *b_mu_msOuterMatchChi2;   //!
   TBranch        *b_mu_msInnerMatchDOF;   //!
   TBranch        *b_mu_msOuterMatchDOF;   //!
   TBranch        *b_mu_nOutliersOnTrack;   //!
   TBranch        *b_mu_nBLHits;   //!
   TBranch        *b_mu_nPixHits;   //!
   TBranch        *b_mu_nSCTHits;   //!
   TBranch        *b_mu_nTRTHits;   //!
   TBranch        *b_mu_nTRTHighTHits;   //!
   TBranch        *b_mu_nBLSharedHits;   //!
   TBranch        *b_mu_nPixSharedHits;   //!
   TBranch        *b_mu_nPixHoles;   //!
   TBranch        *b_mu_nSCTSharedHits;   //!
   TBranch        *b_mu_nSCTHoles;   //!
   TBranch        *b_mu_nTRTOutliers;   //!
   TBranch        *b_mu_nTRTHighTOutliers;   //!
   TBranch        *b_mu_nGangedPixels;   //!
   TBranch        *b_mu_nPixelDeadSensors;   //!
   TBranch        *b_mu_nSCTDeadSensors;   //!
   TBranch        *b_mu_nTRTDeadStraws;   //!
   TBranch        *b_mu_expectBLayerHit;   //!
   TBranch        *b_mu_nPrecisionLayers;   //!
   TBranch        *b_mu_nPrecisionHoleLayers;   //!
   TBranch        *b_mu_nPhiLayers;   //!
   TBranch        *b_mu_nPhiHoleLayers;   //!
   TBranch        *b_mu_nTrigEtaLayers;   //!
   TBranch        *b_mu_nTrigEtaHoleLayers;   //!
   TBranch        *b_mu_primarySector;   //!
   TBranch        *b_mu_secondarySector;   //!
   TBranch        *b_mu_nInnerSmallHits;   //!
   TBranch        *b_mu_nInnerLargeHits;   //!
   TBranch        *b_mu_nMiddleSmallHits;   //!
   TBranch        *b_mu_nMiddleLargeHits;   //!
   TBranch        *b_mu_nOuterSmallHits;   //!
   TBranch        *b_mu_nOuterLargeHits;   //!
   TBranch        *b_mu_nExtendedSmallHits;   //!
   TBranch        *b_mu_nExtendedLargeHits;   //!
   TBranch        *b_mu_nInnerSmallHoles;   //!
   TBranch        *b_mu_nInnerLargeHoles;   //!
   TBranch        *b_mu_nMiddleSmallHoles;   //!
   TBranch        *b_mu_nMiddleLargeHoles;   //!
   TBranch        *b_mu_nOuterSmallHoles;   //!
   TBranch        *b_mu_nOuterLargeHoles;   //!
   TBranch        *b_mu_nExtendedSmallHoles;   //!
   TBranch        *b_mu_nExtendedLargeHoles;   //!
   TBranch        *b_mu_nPhiLayer1Hits;   //!
   TBranch        *b_mu_nPhiLayer2Hits;   //!
   TBranch        *b_mu_nPhiLayer3Hits;   //!
   TBranch        *b_mu_nPhiLayer4Hits;   //!
   TBranch        *b_mu_nEtaLayer1Hits;   //!
   TBranch        *b_mu_nEtaLayer2Hits;   //!
   TBranch        *b_mu_nEtaLayer3Hits;   //!
   TBranch        *b_mu_nEtaLayer4Hits;   //!
   TBranch        *b_mu_nPhiLayer1Holes;   //!
   TBranch        *b_mu_nPhiLayer2Holes;   //!
   TBranch        *b_mu_nPhiLayer3Holes;   //!
   TBranch        *b_mu_nPhiLayer4Holes;   //!
   TBranch        *b_mu_nEtaLayer1Holes;   //!
   TBranch        *b_mu_nEtaLayer2Holes;   //!
   TBranch        *b_mu_nEtaLayer3Holes;   //!
   TBranch        *b_mu_nEtaLayer4Holes;   //!
   TBranch        *b_mu_cb_d0;   //!
   TBranch        *b_mu_cb_z0;   //!
   TBranch        *b_mu_cb_phi0;   //!
   TBranch        *b_mu_cb_theta;   //!
   TBranch        *b_mu_cb_qOverP;   //!
   TBranch        *b_mu_cb_vx;   //!
   TBranch        *b_mu_cb_vy;   //!
   TBranch        *b_mu_cb_vz;   //!
   TBranch        *b_TGC_museg_n;   //!
   TBranch        *b_museg_x;   //!
   TBranch        *b_museg_y;   //!
   TBranch        *b_museg_z;   //!
   TBranch        *b_museg_px;   //!
   TBranch        *b_museg_py;   //!
   TBranch        *b_museg_pz;   //!
   TBranch        *b_museg_t0;   //!
   TBranch        *b_museg_t0error;   //!
   TBranch        *b_museg_chi2;   //!
   TBranch        *b_museg_ndof;   //!
   TBranch        *b_museg_sector;   //!
   TBranch        *b_museg_stationName;   //!
   TBranch        *b_museg_stationEta;   //!
   TBranch        *b_museg_author;   //!
   TBranch        *b_ext_mu_bias_n;   //!
   TBranch        *b_ext_mu_bias_type;   //!
   TBranch        *b_ext_mu_bias_index;   //!
   TBranch        *b_ext_mu_bias_size;   //!
   TBranch        *b_ext_mu_bias_targetVec;   //!
   TBranch        *b_ext_mu_bias_targetDistanceVec;   //!
   TBranch        *b_ext_mu_bias_targetEtaVec;   //!
   TBranch        *b_ext_mu_bias_targetPhiVec;   //!
   TBranch        *b_ext_mu_bias_targetDeltaEtaVec;   //!
   TBranch        *b_ext_mu_bias_targetDeltaPhiVec;   //!
   TBranch        *b_ext_mu_bias_targetPxVec;   //!
   TBranch        *b_ext_mu_bias_targetPyVec;   //!
   TBranch        *b_ext_mu_bias_targetPzVec;   //!
   TBranch        *b_ext_mu_ubias_n;   //!
   TBranch        *b_ext_mu_ubias_type;   //!
   TBranch        *b_ext_mu_ubias_index;   //!
   TBranch        *b_ext_mu_ubias_size;   //!
   TBranch        *b_ext_mu_ubias_targetVec;   //!
   TBranch        *b_ext_mu_ubias_targetDistanceVec;   //!
   TBranch        *b_ext_mu_ubias_targetEtaVec;   //!
   TBranch        *b_ext_mu_ubias_targetPhiVec;   //!
   TBranch        *b_ext_mu_ubias_targetDeltaEtaVec;   //!
   TBranch        *b_ext_mu_ubias_targetDeltaPhiVec;   //!
   TBranch        *b_ext_mu_ubias_targetPxVec;   //!
   TBranch        *b_ext_mu_ubias_targetPyVec;   //!
   TBranch        *b_ext_mu_ubias_targetPzVec;   //!
   TBranch        *b_trigger_info_n;   //!
   TBranch        *b_trigger_info_chain;   //!
   TBranch        *b_trigger_info_isPassed;   //!
   TBranch        *b_trigger_info_nTracks;   //!
   TBranch        *b_trigger_info_typeVec;   //!
   TBranch        *b_trigger_info_ptVec;   //!
   TBranch        *b_trigger_info_etaVec;   //!
   TBranch        *b_trigger_info_phiVec;   //!
   TBranch        *b_vxp_n;   //!
   TBranch        *b_vxp_x;   //!
   TBranch        *b_vxp_y;   //!
   TBranch        *b_vxp_z;   //!
   TBranch        *b_vxp_cov_x;   //!
   TBranch        *b_vxp_cov_y;   //!
   TBranch        *b_vxp_cov_z;   //!
   TBranch        *b_vxp_cov_xy;   //!
   TBranch        *b_vxp_cov_xz;   //!
   TBranch        *b_vxp_cov_yz;   //!
   TBranch        *b_vxp_chi2;   //!
   TBranch        *b_vxp_ndof;   //!
   TBranch        *b_vxp_nTracks;   //!
   TBranch        *b_vxp_type;   //!
   TBranch        *b_TGC_prd_n;   //!
   TBranch        *b_TGC_prd_x;   //!
   TBranch        *b_TGC_prd_y;   //!
   TBranch        *b_TGC_prd_z;   //!
   TBranch        *b_TGC_prd_shortWidth;   //!
   TBranch        *b_TGC_prd_longWidth;   //!
   TBranch        *b_TGC_prd_length;   //!
   TBranch        *b_TGC_prd_isStrip;   //!
   TBranch        *b_TGC_prd_gasGap;   //!
   TBranch        *b_TGC_prd_channel;   //!
   TBranch        *b_TGC_prd_eta;   //!
   TBranch        *b_TGC_prd_phi;   //!
   TBranch        *b_TGC_prd_station;   //!
   TBranch        *b_TGC_prd_bunch;   //!
   TBranch        *b_MDT_prd_n;   //!
   TBranch        *b_MDT_prd_x;   //!
   TBranch        *b_MDT_prd_y;   //!
   TBranch        *b_MDT_prd_z;   //!
   TBranch        *b_MDT_prd_adc;   //!
   TBranch        *b_MDT_prd_tdc;   //!
   TBranch        *b_MDT_prd_status;   //!
   TBranch        *b_MDT_prd_drift_radius;   //!
   TBranch        *b_MDT_prd_drift_radius_error;   //!
   TBranch        *b_TGC_coin_n;   //!
   TBranch        *b_TGC_coin_x_In;   //!
   TBranch        *b_TGC_coin_y_In;   //!
   TBranch        *b_TGC_coin_z_In;   //!
   TBranch        *b_TGC_coin_x_Out;   //!
   TBranch        *b_TGC_coin_y_Out;   //!
   TBranch        *b_TGC_coin_z_Out;   //!
   TBranch        *b_TGC_coin_width_In;   //!
   TBranch        *b_TGC_coin_width_Out;   //!
   TBranch        *b_TGC_coin_width_R;   //!
   TBranch        *b_TGC_coin_width_Phi;   //!
   TBranch        *b_TGC_coin_isAside;   //!
   TBranch        *b_TGC_coin_isForward;   //!
   TBranch        *b_TGC_coin_isStrip;   //!
   TBranch        *b_TGC_coin_isInner;   //!
   TBranch        *b_TGC_coin_isPositiveDeltaR;   //!
   TBranch        *b_TGC_coin_type;   //!
   TBranch        *b_TGC_coin_trackletId;   //!
   TBranch        *b_TGC_coin_trackletIdStrip;   //!
   TBranch        *b_TGC_coin_phi;   //!
   TBranch        *b_TGC_coin_roi;   //!
   TBranch        *b_TGC_coin_pt;   //!
   TBranch        *b_TGC_coin_delta;   //!
   TBranch        *b_TGC_coin_sub;   //!
   TBranch        *b_TGC_coin_veto;   //!
   TBranch        *b_TGC_coin_bunch;   //!
   TBranch        *b_TGC_coin_inner;   //!
   TBranch        *b_TILE_murcv_trig_n;   //!
   TBranch        *b_TILE_murcv_trig_mod;   //!
   TBranch        *b_TILE_murcv_trig_part;   //!
   TBranch        *b_TILE_murcv_trig_bit0;   //!
   TBranch        *b_TILE_murcv_trig_bit1;   //!
   TBranch        *b_TILE_murcv_trig_bit2;   //!
   TBranch        *b_TILE_murcv_trig_bit3;   //!
   TBranch        *b_TILE_murcv_raw_n;   //!
   TBranch        *b_TILE_murcv_raw_count;   //!
   TBranch        *b_TILE_murcv_raw_energy;   //!
   TBranch        *b_TILE_murcv_raw_ros;   //!
   TBranch        *b_TILE_murcv_raw_drawer;   //!
   TBranch        *b_TILE_murcv_raw_channel;   //!
   TBranch        *b_TILE_murcv_digit_n;   //!
   TBranch        *b_TILE_murcv_digit_nSamples;   //!
   TBranch        *b_TILE_murcv_digit_ros;   //!
   TBranch        *b_TILE_murcv_digit_drawer;   //!
   TBranch        *b_TILE_murcv_digit_channel;   //!
   TBranch        *b_TILE_murcv_digit_sampleVec;   //!
   TBranch        *b_TGC_hierarchy_n;   //!
   TBranch        *b_TGC_hierarchy_index;   //!
   TBranch        *b_TGC_hierarchy_dR_hiPt;   //!
   TBranch        *b_TGC_hierarchy_dPhi_hiPt;   //!
   TBranch        *b_TGC_hierarchy_dR_tracklet;   //!
   TBranch        *b_TGC_hierarchy_dPhi_tracklet;   //!
   TBranch        *b_TGC_hierarchy_isChamberBoundary;   //!
   TBranch        *b_muctpi_candidateMultiplicities;   //!
   TBranch        *b_muctpi_nDataWords;   //!
   TBranch        *b_muctpi_dataWords;   //!
   TBranch        *b_muctpi_dw_eta;   //!
   TBranch        *b_muctpi_dw_phi;   //!
   TBranch        *b_muctpi_dw_source;   //!
   TBranch        *b_muctpi_dw_hemisphere;   //!
   TBranch        *b_muctpi_dw_bcid;   //!
   TBranch        *b_muctpi_dw_sectorID;   //!
   TBranch        *b_muctpi_dw_thrNumber;   //!
   TBranch        *b_muctpi_dw_roi;   //!
   TBranch        *b_muctpi_dw_veto;   //!
   TBranch        *b_muctpi_dw_firstCandidate;   //!
   TBranch        *b_muctpi_dw_moreCandInRoI;   //!
   TBranch        *b_muctpi_dw_moreCandInSector;   //!
   TBranch        *b_muctpi_dw_charge;   //!
   TBranch        *b_muctpi_dw_candidateVetoed;   //!
  //  TBranch        *b_TILE_cell_n;   //!
  //  TBranch        *b_TILE_cell_E;   //!
  //  TBranch        *b_TILE_cell_Et;   //!
  //  TBranch        *b_TILE_cell_eta;   //!
  //  TBranch        *b_TILE_cell_phi;   //!
  //  TBranch        *b_TILE_cell_sinTh;   //!
  //  TBranch        *b_TILE_cell_cosTh;   //!
  //  TBranch        *b_TILE_cell_cotTh;   //!
  //  TBranch        *b_TILE_cell_x;   //!
  //  TBranch        *b_TILE_cell_y;   //!
  //  TBranch        *b_TILE_cell_z;   //!
  //  TBranch        *b_TILE_cell_badcell;   //!
  //  TBranch        *b_TILE_cell_partition;   //!
  //  TBranch        *b_TILE_cell_section;   //!
  //  TBranch        *b_TILE_cell_side;   //!
  //  TBranch        *b_TILE_cell_module;   //!
  //  TBranch        *b_TILE_cell_tower;   //!
  //  TBranch        *b_TILE_cell_sample;   //!

   TMDB(TTree *tree1=0);
   virtual ~TMDB();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntryPhysics(Long64_t entry);
   virtual Long64_t LoadTreePhysics(Long64_t entry);
   virtual void     Init(TTree *tree1);
   virtual void     LoopPhysics(string fout);
   virtual Bool_t   Notify();
   virtual void     ShowPhysics(Long64_t entry = -1);

//   virtual void tgc_trigger_loop(bool ignoreAside);
//   virtual void coin_loop();
//   virtual void raw_loop();
//   virtual void trig_loop();
//   virtual void tmdb_sent();
//   virtual void tmdb_digit_loop(Int_t module_num, string filename, bool ignoreAside);
//   virtual void sector_loop();
//   virtual void module_loop();
//   virtual void tile_offline_loop(bool ignoreAside);
//   virtual void fit_loop();

//   virtual Float_t calc_dR(Float_t dEta, Float_t dPhi);

//   virtual void  UserInit();
//   virtual void  UserFin(string fout);
//   virtual Int_t tilemodule_to_check(Int_t sector);
//   virtual Int_t modulenumber_in_glink(Int_t sector);

	virtual void 		UserInit();
	virtual void		loop_offline(bool ignoreAside);
  virtual void		loop_hlt(bool ignoreAside);
  virtual void		loop_l1_trigger(bool ignoreAside);
  virtual void		loop_tmdb(bool ignoreAside);
  virtual void		loop_extra(bool ignoreAside);
	virtual Int_t		tile_phi_module(Float_t phi);
	virtual Float_t 	calc_dR(Float_t dEta, Float_t dPhi);

  virtual void		print_root(string fout);
  virtual Int_t   tilemodule_to_check(Int_t sector);


  // Variables
	Int_t 					npulses[2][total_tile_modules][total_tmdb_mod_ch];

  TH1**           th1_offline;
  TH1**           th1_offline_eta;
  TH1**           th1_offline_phi;
  TH1**           th1_offline_pt;

  // TH1**           th1_hlt;
  TH1**           th1_hlt_eta;
  TH1**           th1_hlt_eta_vp;
  TH1**           th1_hlt_eta_vn;
  TH1**           th1_hlt_eta_fp;
  TH1**           th1_hlt_eta_fn;
  TH1**           th1_hlt_phi;
  TH1**           th1_hlt_phi_vp;
  TH1**           th1_hlt_phi_vn;
  TH1**           th1_hlt_phi_fp;
  TH1**           th1_hlt_phi_fn;
  TH1**           th1_hlt_pt;
  TH1**           th1_hlt_pt_vp;
  TH1**           th1_hlt_pt_vn;
  TH1**           th1_hlt_pt_fp;
  TH1**           th1_hlt_pt_fn;

  TH1*            th1_l1muon;
  TH1*            th1_l1muon_eta;
  TH1*            th1_l1muon_phi;
  // TH1*            th1_l1muon_pt;

  // TH1*            th1_l1muon_15;
  TH1*            th1_l1muon_15_eta;
  TH1**           th1_l1muon_15_eta_vp;
  TH1**           th1_l1muon_15_eta_vn;
  TH1**           th1_l1muon_15_eta_fp;
  TH1**           th1_l1muon_15_eta_fn;
  TH1**           th1_l1muon_15_phi;
  TH1**           th1_l1muon_15_phi_vp;
  TH1**           th1_l1muon_15_phi_vn;
  TH1**           th1_l1muon_15_phi_fp;
  TH1**           th1_l1muon_15_phi_fn;
  TH1**           th1_l1muon_15_pt;
  TH1**           th1_l1muon_15_pt_vp;
  TH1**           th1_l1muon_15_pt_vn;
  TH1**           th1_l1muon_15_pt_fp;
  TH1**           th1_l1muon_15_pt_fn;
  // TH1*            th1_l1muon_20;
  TH1*            th1_l1muon_20_eta;
  TH1**           th1_l1muon_20_eta_vp;
  TH1**           th1_l1muon_20_eta_vn;
  TH1**           th1_l1muon_20_eta_fp;
  TH1**           th1_l1muon_20_eta_fn;
  TH1**           th1_l1muon_20_phi;
  TH1**           th1_l1muon_20_phi_vp;
  TH1**           th1_l1muon_20_phi_vn;
  TH1**           th1_l1muon_20_phi_fp;
  TH1**           th1_l1muon_20_phi_fn;
  TH1**           th1_l1muon_20_pt;
  TH1**           th1_l1muon_20_pt_vp;
  TH1**           th1_l1muon_20_pt_vn;
  TH1**           th1_l1muon_20_pt_fp;
  TH1**           th1_l1muon_20_pt_fn;

  TH1*            th1_tmdb;
  TH1*            th1_tmdb_eta;
  TH1**           th1_tmdb_phi;
  // TH1*            th1_tmdb_15;
  // TH1*            th1_tmdb_15_eta;
  TH1**           th1_tmdb_15_eta_vp;
  TH1**           th1_tmdb_15_eta_vn;
  TH1**           th1_tmdb_15_eta_fp;
  TH1**           th1_tmdb_15_eta_fn;
  // TH1**           th1_tmdb_15_phi;
  TH1**           th1_tmdb_15_phi_vp;
  TH1**           th1_tmdb_15_phi_vn;
  TH1**           th1_tmdb_15_phi_fp;
  TH1**           th1_tmdb_15_phi_fn;
  TH1**           th1_tmdb_15_pt;
  TH1**           th1_tmdb_15_pt_vp;
  TH1**           th1_tmdb_15_pt_vn;
  TH1**           th1_tmdb_15_pt_fp;
  TH1**           th1_tmdb_15_pt_fn;
  // TH1*            th1_tmdb_20;
  // TH1*            th1_tmdb_20_eta;
  TH1**           th1_tmdb_20_eta_vp;
  TH1**           th1_tmdb_20_eta_vn;
  TH1**           th1_tmdb_20_eta_fp;
  TH1**           th1_tmdb_20_eta_fn;
  // TH1**           th1_tmdb_20_phi;
  TH1**           th1_tmdb_20_phi_vp;
  TH1**           th1_tmdb_20_phi_vn;
  TH1**           th1_tmdb_20_phi_fp;
  TH1**           th1_tmdb_20_phi_fn;
  TH1**           th1_tmdb_20_pt;
  TH1**           th1_tmdb_20_pt_vp;
  TH1**           th1_tmdb_20_pt_vn;
  TH1**           th1_tmdb_20_pt_fp;
  TH1**           th1_tmdb_20_pt_fn;

  TH1**           th1_mf_out;
  TH1**           th1_mf_decision;
  TH1**           th1_mf_energy;

};

#endif

#ifdef TMDB_cxx
TMDB::TMDB(TTree *tree1) : fChainPhysics(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
	// Tree from L1TGCNtuple
   if (tree1 == 0) {
      TFile *f1 = (TFile*)gROOT->GetListOfFiles()->FindObject("L1TGCNtuple.root");
      if (!f1 || !f1->IsOpen()) {
         f1 = new TFile("L1TGCNtuple.root");
      }
      f1->GetObject("physics",tree1);

   }

   Init(tree1);
}

TMDB::~TMDB()
{
   if (!fChainPhysics) return;
   delete fChainPhysics->GetCurrentFile();
}

Int_t TMDB::GetEntryPhysics(Long64_t entry)
{
// Read contents of L1TGCNtuple entry.
   if (!fChainPhysics) return 0;
   return fChainPhysics->GetEntry(entry);
}

Long64_t TMDB::LoadTreePhysics(Long64_t entry)
{
// Set the environment to read one L1TGCNtuple entry
   if (!fChainPhysics) return -5;
   Long64_t centry = fChainPhysics->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChainPhysics->GetTreeNumber() != fCurrentPhysics) {
      fCurrentPhysics = fChainPhysics->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TMDB::Init(TTree *tree1)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer for h2000

   // Set object pointer for L1TGCNtuple
   trig_L1_mu_eta = 0;
   trig_L1_mu_phi = 0;
   trig_L1_mu_thrName = 0;
   trig_L1_mu_thrNumber = 0;
   trig_L1_mu_RoINumber = 0;
   trig_L1_mu_sectorAddress = 0;
   trig_L1_mu_firstCandidate = 0;
   trig_L1_mu_moreCandInRoI = 0;
   trig_L1_mu_moreCandInSector = 0;
   trig_L1_mu_source = 0;
   trig_L1_mu_hemisphere = 0;
   trig_L1_mu_charge = 0;
   trig_L1_mu_vetoed = 0;
   mu_pt = 0;
   mu_eta = 0;
   mu_phi = 0;
   mu_m = 0;
   mu_charge = 0;
   mu_author = 0;
   mu_allAuthors = 0;
   mu_muonType = 0;
   mu_etcone20 = 0;
   mu_etcone30 = 0;
   mu_etcone40 = 0;
   mu_ptcone20 = 0;
   mu_ptcone30 = 0;
   mu_ptcone40 = 0;
   mu_trackfitchi2 = 0;
   mu_trackfitndof = 0;
   mu_msInnerMatchChi2 = 0;
   mu_msOuterMatchChi2 = 0;
   mu_msInnerMatchDOF = 0;
   mu_msOuterMatchDOF = 0;
   mu_nOutliersOnTrack = 0;
   mu_nBLHits = 0;
   mu_nPixHits = 0;
   mu_nSCTHits = 0;
   mu_nTRTHits = 0;
   mu_nTRTHighTHits = 0;
   mu_nBLSharedHits = 0;
   mu_nPixSharedHits = 0;
   mu_nPixHoles = 0;
   mu_nSCTSharedHits = 0;
   mu_nSCTHoles = 0;
   mu_nTRTOutliers = 0;
   mu_nTRTHighTOutliers = 0;
   mu_nGangedPixels = 0;
   mu_nPixelDeadSensors = 0;
   mu_nSCTDeadSensors = 0;
   mu_nTRTDeadStraws = 0;
   mu_expectBLayerHit = 0;
   mu_nPrecisionLayers = 0;
   mu_nPrecisionHoleLayers = 0;
   mu_nPhiLayers = 0;
   mu_nPhiHoleLayers = 0;
   mu_nTrigEtaLayers = 0;
   mu_nTrigEtaHoleLayers = 0;
   mu_primarySector = 0;
   mu_secondarySector = 0;
   mu_nInnerSmallHits = 0;
   mu_nInnerLargeHits = 0;
   mu_nMiddleSmallHits = 0;
   mu_nMiddleLargeHits = 0;
   mu_nOuterSmallHits = 0;
   mu_nOuterLargeHits = 0;
   mu_nExtendedSmallHits = 0;
   mu_nExtendedLargeHits = 0;
   mu_nInnerSmallHoles = 0;
   mu_nInnerLargeHoles = 0;
   mu_nMiddleSmallHoles = 0;
   mu_nMiddleLargeHoles = 0;
   mu_nOuterSmallHoles = 0;
   mu_nOuterLargeHoles = 0;
   mu_nExtendedSmallHoles = 0;
   mu_nExtendedLargeHoles = 0;
   mu_nPhiLayer1Hits = 0;
   mu_nPhiLayer2Hits = 0;
   mu_nPhiLayer3Hits = 0;
   mu_nPhiLayer4Hits = 0;
   mu_nEtaLayer1Hits = 0;
   mu_nEtaLayer2Hits = 0;
   mu_nEtaLayer3Hits = 0;
   mu_nEtaLayer4Hits = 0;
   mu_nPhiLayer1Holes = 0;
   mu_nPhiLayer2Holes = 0;
   mu_nPhiLayer3Holes = 0;
   mu_nPhiLayer4Holes = 0;
   mu_nEtaLayer1Holes = 0;
   mu_nEtaLayer2Holes = 0;
   mu_nEtaLayer3Holes = 0;
   mu_nEtaLayer4Holes = 0;
   mu_cb_d0 = 0;
   mu_cb_z0 = 0;
   mu_cb_phi0 = 0;
   mu_cb_theta = 0;
   mu_cb_qOverP = 0;
   mu_cb_vx = 0;
   mu_cb_vy = 0;
   mu_cb_vz = 0;
   museg_x = 0;
   museg_y = 0;
   museg_z = 0;
   museg_px = 0;
   museg_py = 0;
   museg_pz = 0;
   museg_t0 = 0;
   museg_t0error = 0;
   museg_chi2 = 0;
   museg_ndof = 0;
   museg_sector = 0;
   museg_stationName = 0;
   museg_stationEta = 0;
   museg_author = 0;
   ext_mu_bias_type = 0;
   ext_mu_bias_index = 0;
   ext_mu_bias_size = 0;
   ext_mu_bias_targetVec = 0;
   ext_mu_bias_targetDistanceVec = 0;
   ext_mu_bias_targetEtaVec = 0;
   ext_mu_bias_targetPhiVec = 0;
   ext_mu_bias_targetDeltaEtaVec = 0;
   ext_mu_bias_targetDeltaPhiVec = 0;
   ext_mu_bias_targetPxVec = 0;
   ext_mu_bias_targetPyVec = 0;
   ext_mu_bias_targetPzVec = 0;
   ext_mu_ubias_type = 0;
   ext_mu_ubias_index = 0;
   ext_mu_ubias_size = 0;
   ext_mu_ubias_targetVec = 0;
   ext_mu_ubias_targetDistanceVec = 0;
   ext_mu_ubias_targetEtaVec = 0;
   ext_mu_ubias_targetPhiVec = 0;
   ext_mu_ubias_targetDeltaEtaVec = 0;
   ext_mu_ubias_targetDeltaPhiVec = 0;
   ext_mu_ubias_targetPxVec = 0;
   ext_mu_ubias_targetPyVec = 0;
   ext_mu_ubias_targetPzVec = 0;
   trigger_info_chain = 0;
   trigger_info_isPassed = 0;
   trigger_info_nTracks = 0;
   trigger_info_typeVec = 0;
   trigger_info_ptVec = 0;
   trigger_info_etaVec = 0;
   trigger_info_phiVec = 0;
   vxp_x = 0;
   vxp_y = 0;
   vxp_z = 0;
   vxp_cov_x = 0;
   vxp_cov_y = 0;
   vxp_cov_z = 0;
   vxp_cov_xy = 0;
   vxp_cov_xz = 0;
   vxp_cov_yz = 0;
   vxp_chi2 = 0;
   vxp_ndof = 0;
   vxp_nTracks = 0;
   vxp_type = 0;
   TGC_prd_x = 0;
   TGC_prd_y = 0;
   TGC_prd_z = 0;
   TGC_prd_shortWidth = 0;
   TGC_prd_longWidth = 0;
   TGC_prd_length = 0;
   TGC_prd_isStrip = 0;
   TGC_prd_gasGap = 0;
   TGC_prd_channel = 0;
   TGC_prd_eta = 0;
   TGC_prd_phi = 0;
   TGC_prd_station = 0;
   TGC_prd_bunch = 0;
   MDT_prd_x = 0;
   MDT_prd_y = 0;
   MDT_prd_z = 0;
   MDT_prd_adc = 0;
   MDT_prd_tdc = 0;
   MDT_prd_status = 0;
   MDT_prd_drift_radius = 0;
   MDT_prd_drift_radius_error = 0;
   TGC_coin_x_In = 0;
   TGC_coin_y_In = 0;
   TGC_coin_z_In = 0;
   TGC_coin_x_Out = 0;
   TGC_coin_y_Out = 0;
   TGC_coin_z_Out = 0;
   TGC_coin_width_In = 0;
   TGC_coin_width_Out = 0;
   TGC_coin_width_R = 0;
   TGC_coin_width_Phi = 0;
   TGC_coin_isAside = 0;
   TGC_coin_isForward = 0;
   TGC_coin_isStrip = 0;
   TGC_coin_isInner = 0;
   TGC_coin_isPositiveDeltaR = 0;
   TGC_coin_type = 0;
   TGC_coin_trackletId = 0;
   TGC_coin_trackletIdStrip = 0;
   TGC_coin_phi = 0;
   TGC_coin_roi = 0;
   TGC_coin_pt = 0;
   TGC_coin_delta = 0;
   TGC_coin_sub = 0;
   TGC_coin_veto = 0;
   TGC_coin_bunch = 0;
   TGC_coin_inner = 0;
   TILE_murcv_trig_mod = 0;
   TILE_murcv_trig_part = 0;
   TILE_murcv_trig_bit0 = 0;
   TILE_murcv_trig_bit1 = 0;
   TILE_murcv_trig_bit2 = 0;
   TILE_murcv_trig_bit3 = 0;
   TILE_murcv_raw_count = 0;
   TILE_murcv_raw_energy = 0;
   TILE_murcv_raw_ros = 0;
   TILE_murcv_raw_drawer = 0;
   TILE_murcv_raw_channel = 0;
   TILE_murcv_digit_nSamples = 0;
   TILE_murcv_digit_ros = 0;
   TILE_murcv_digit_drawer = 0;
   TILE_murcv_digit_channel = 0;
   TILE_murcv_digit_sampleVec = 0;
   TGC_hierarchy_index = 0;
   TGC_hierarchy_dR_hiPt = 0;
   TGC_hierarchy_dPhi_hiPt = 0;
   TGC_hierarchy_dR_tracklet = 0;
   TGC_hierarchy_dPhi_tracklet = 0;
   TGC_hierarchy_isChamberBoundary = 0;
   muctpi_candidateMultiplicities = 0;
   muctpi_dataWords = 0;
   muctpi_dw_eta = 0;
   muctpi_dw_phi = 0;
   muctpi_dw_source = 0;
   muctpi_dw_hemisphere = 0;
   muctpi_dw_bcid = 0;
   muctpi_dw_sectorID = 0;
   muctpi_dw_thrNumber = 0;
   muctpi_dw_roi = 0;
   muctpi_dw_veto = 0;
   muctpi_dw_firstCandidate = 0;
   muctpi_dw_moreCandInRoI = 0;
   muctpi_dw_moreCandInSector = 0;
   muctpi_dw_charge = 0;
   muctpi_dw_candidateVetoed = 0;
  //  TILE_cell_E = 0;
  //  TILE_cell_Et = 0;
  //  TILE_cell_eta = 0;
  //  TILE_cell_phi = 0;
  //  TILE_cell_sinTh = 0;
  //  TILE_cell_cosTh = 0;
  //  TILE_cell_cotTh = 0;
  //  TILE_cell_x = 0;
  //  TILE_cell_y = 0;
  //  TILE_cell_z = 0;
  //  TILE_cell_badcell = 0;
  //  TILE_cell_partition = 0;
  //  TILE_cell_section = 0;
  //  TILE_cell_side = 0;
  //  TILE_cell_module = 0;
  //  TILE_cell_tower = 0;
  //  TILE_cell_sample = 0;
   // Set branch addresses and branch pointers for L1TGCNtuple
   if (!tree1) return;
   fChainPhysics = tree1;
   fCurrentPhysics = -1;
   fChainPhysics->SetMakeClass(1);

   fChainPhysics->SetBranchAddress("RunNumber", &RunNumber, &b_runNumber);
   fChainPhysics->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChainPhysics->SetBranchAddress("timeStamp", &timeStamp, &b_timeStamp);
   fChainPhysics->SetBranchAddress("timeStampNSOffset", &timeStampNSOffset, &b_timeStampNSOffset);
   fChainPhysics->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChainPhysics->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChainPhysics->SetBranchAddress("detmask0", &detmask0, &b_detmask0);
   fChainPhysics->SetBranchAddress("detmask1", &detmask1, &b_detmask1);
   fChainPhysics->SetBranchAddress("actualIntPerXing", &actualIntPerXing, &b_actualIntPerXing);
   fChainPhysics->SetBranchAddress("averageIntPerXing", &averageIntPerXing, &b_averageIntPerXing);
   fChainPhysics->SetBranchAddress("pixelFlags", &pixelFlags, &b_pixelFlags);
   fChainPhysics->SetBranchAddress("sctFlags", &sctFlags, &b_sctFlags);
   fChainPhysics->SetBranchAddress("trtFlags", &trtFlags, &b_trtFlags);
   fChainPhysics->SetBranchAddress("larFlags", &larFlags, &b_larFlags);
   fChainPhysics->SetBranchAddress("tileFlags", &tileFlags, &b_tileFlags);
   fChainPhysics->SetBranchAddress("muonFlags", &muonFlags, &b_muonFlags);
   fChainPhysics->SetBranchAddress("fwdFlags", &fwdFlags, &b_fwdFlags);
   fChainPhysics->SetBranchAddress("coreFlags", &coreFlags, &b_coreFlags);
   fChainPhysics->SetBranchAddress("pixelError", &pixelError, &b_pixelError);
   fChainPhysics->SetBranchAddress("sctError", &sctError, &b_sctError);
   fChainPhysics->SetBranchAddress("trtError", &trtError, &b_trtError);
   fChainPhysics->SetBranchAddress("larError", &larError, &b_larError);
   fChainPhysics->SetBranchAddress("tileError", &tileError, &b_tileError);
   fChainPhysics->SetBranchAddress("muonError", &muonError, &b_muonError);
   fChainPhysics->SetBranchAddress("fwdError", &fwdError, &b_fwdError);
   fChainPhysics->SetBranchAddress("coreError", &coreError, &b_coreError);
   fChainPhysics->SetBranchAddress("trig_L1_mu_n", &trig_L1_mu_n, &b_trig_L1_mu_n);
   fChainPhysics->SetBranchAddress("trig_L1_mu_eta", &trig_L1_mu_eta, &b_trig_L1_mu_eta);
   fChainPhysics->SetBranchAddress("trig_L1_mu_phi", &trig_L1_mu_phi, &b_trig_L1_mu_phi);
   fChainPhysics->SetBranchAddress("trig_L1_mu_thrName", &trig_L1_mu_thrName, &b_trig_L1_mu_thrName);
   fChainPhysics->SetBranchAddress("trig_L1_mu_thrNumber", &trig_L1_mu_thrNumber, &b_trig_L1_mu_thrNumber);
   fChainPhysics->SetBranchAddress("trig_L1_mu_RoINumber", &trig_L1_mu_RoINumber, &b_trig_L1_mu_RoINumber);
   fChainPhysics->SetBranchAddress("trig_L1_mu_sectorAddress", &trig_L1_mu_sectorAddress, &b_trig_L1_mu_sectorAddress);
   fChainPhysics->SetBranchAddress("trig_L1_mu_firstCandidate", &trig_L1_mu_firstCandidate, &b_trig_L1_mu_firstCandidate);
   fChainPhysics->SetBranchAddress("trig_L1_mu_moreCandInRoI", &trig_L1_mu_moreCandInRoI, &b_trig_L1_mu_moreCandInRoI);
   fChainPhysics->SetBranchAddress("trig_L1_mu_moreCandInSector", &trig_L1_mu_moreCandInSector, &b_trig_L1_mu_moreCandInSector);
   fChainPhysics->SetBranchAddress("trig_L1_mu_source", &trig_L1_mu_source, &b_trig_L1_mu_source);
   fChainPhysics->SetBranchAddress("trig_L1_mu_hemisphere", &trig_L1_mu_hemisphere, &b_trig_L1_mu_hemisphere);
   fChainPhysics->SetBranchAddress("trig_L1_mu_charge", &trig_L1_mu_charge, &b_trig_L1_mu_charge);
   fChainPhysics->SetBranchAddress("trig_L1_mu_vetoed", &trig_L1_mu_vetoed, &b_trig_L1_mu_vetoed);
   fChainPhysics->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
   fChainPhysics->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   fChainPhysics->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChainPhysics->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChainPhysics->SetBranchAddress("mu_m", &mu_m, &b_mu_m);
   fChainPhysics->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChainPhysics->SetBranchAddress("mu_author", &mu_author, &b_mu_author);
   fChainPhysics->SetBranchAddress("mu_allAuthors", &mu_allAuthors, &b_mu_allAuthors);
   fChainPhysics->SetBranchAddress("mu_muonType", &mu_muonType, &b_mu_muonType);
   fChainPhysics->SetBranchAddress("mu_etcone20", &mu_etcone20, &b_mu_etcone20);
   fChainPhysics->SetBranchAddress("mu_etcone30", &mu_etcone30, &b_mu_etcone30);
   fChainPhysics->SetBranchAddress("mu_etcone40", &mu_etcone40, &b_mu_etcone40);
   fChainPhysics->SetBranchAddress("mu_ptcone20", &mu_ptcone20, &b_mu_ptcone20);
   fChainPhysics->SetBranchAddress("mu_ptcone30", &mu_ptcone30, &b_mu_ptcone30);
   fChainPhysics->SetBranchAddress("mu_ptcone40", &mu_ptcone40, &b_mu_ptcone40);
   fChainPhysics->SetBranchAddress("mu_trackfitchi2", &mu_trackfitchi2, &b_mu_trackfitchi2);
   fChainPhysics->SetBranchAddress("mu_trackfitndof", &mu_trackfitndof, &b_mu_trackfitndof);
   fChainPhysics->SetBranchAddress("mu_isPassedMCP", &mu_isPassedMCP, &b_mu_isPassedMCP);
   fChainPhysics->SetBranchAddress("mu_quality", &mu_quality, &b_mu_quality);
   fChainPhysics->SetBranchAddress("mu_msInnerMatchChi2", &mu_msInnerMatchChi2, &b_mu_msInnerMatchChi2);
   fChainPhysics->SetBranchAddress("mu_msOuterMatchChi2", &mu_msOuterMatchChi2, &b_mu_msOuterMatchChi2);
   fChainPhysics->SetBranchAddress("mu_msInnerMatchDOF", &mu_msInnerMatchDOF, &b_mu_msInnerMatchDOF);
   fChainPhysics->SetBranchAddress("mu_msOuterMatchDOF", &mu_msOuterMatchDOF, &b_mu_msOuterMatchDOF);
   fChainPhysics->SetBranchAddress("mu_nOutliersOnTrack", &mu_nOutliersOnTrack, &b_mu_nOutliersOnTrack);
   fChainPhysics->SetBranchAddress("mu_nBLHits", &mu_nBLHits, &b_mu_nBLHits);
   fChainPhysics->SetBranchAddress("mu_nPixHits", &mu_nPixHits, &b_mu_nPixHits);
   fChainPhysics->SetBranchAddress("mu_nSCTHits", &mu_nSCTHits, &b_mu_nSCTHits);
   fChainPhysics->SetBranchAddress("mu_nTRTHits", &mu_nTRTHits, &b_mu_nTRTHits);
   fChainPhysics->SetBranchAddress("mu_nTRTHighTHits", &mu_nTRTHighTHits, &b_mu_nTRTHighTHits);
   fChainPhysics->SetBranchAddress("mu_nBLSharedHits", &mu_nBLSharedHits, &b_mu_nBLSharedHits);
   fChainPhysics->SetBranchAddress("mu_nPixSharedHits", &mu_nPixSharedHits, &b_mu_nPixSharedHits);
   fChainPhysics->SetBranchAddress("mu_nPixHoles", &mu_nPixHoles, &b_mu_nPixHoles);
   fChainPhysics->SetBranchAddress("mu_nSCTSharedHits", &mu_nSCTSharedHits, &b_mu_nSCTSharedHits);
   fChainPhysics->SetBranchAddress("mu_nSCTHoles", &mu_nSCTHoles, &b_mu_nSCTHoles);
   fChainPhysics->SetBranchAddress("mu_nTRTOutliers", &mu_nTRTOutliers, &b_mu_nTRTOutliers);
   fChainPhysics->SetBranchAddress("mu_nTRTHighTOutliers", &mu_nTRTHighTOutliers, &b_mu_nTRTHighTOutliers);
   fChainPhysics->SetBranchAddress("mu_nGangedPixels", &mu_nGangedPixels, &b_mu_nGangedPixels);
   fChainPhysics->SetBranchAddress("mu_nPixelDeadSensors", &mu_nPixelDeadSensors, &b_mu_nPixelDeadSensors);
   fChainPhysics->SetBranchAddress("mu_nSCTDeadSensors", &mu_nSCTDeadSensors, &b_mu_nSCTDeadSensors);
   fChainPhysics->SetBranchAddress("mu_nTRTDeadStraws", &mu_nTRTDeadStraws, &b_mu_nTRTDeadStraws);
   fChainPhysics->SetBranchAddress("mu_expectBLayerHit", &mu_expectBLayerHit, &b_mu_expectBLayerHit);
   fChainPhysics->SetBranchAddress("mu_nPrecisionLayers", &mu_nPrecisionLayers, &b_mu_nPrecisionLayers);
   fChainPhysics->SetBranchAddress("mu_nPrecisionHoleLayers", &mu_nPrecisionHoleLayers, &b_mu_nPrecisionHoleLayers);
   fChainPhysics->SetBranchAddress("mu_nPhiLayers", &mu_nPhiLayers, &b_mu_nPhiLayers);
   fChainPhysics->SetBranchAddress("mu_nPhiHoleLayers", &mu_nPhiHoleLayers, &b_mu_nPhiHoleLayers);
   fChainPhysics->SetBranchAddress("mu_nTrigEtaLayers", &mu_nTrigEtaLayers, &b_mu_nTrigEtaLayers);
   fChainPhysics->SetBranchAddress("mu_nTrigEtaHoleLayers", &mu_nTrigEtaHoleLayers, &b_mu_nTrigEtaHoleLayers);
   fChainPhysics->SetBranchAddress("mu_primarySector", &mu_primarySector, &b_mu_primarySector);
   fChainPhysics->SetBranchAddress("mu_secondarySector", &mu_secondarySector, &b_mu_secondarySector);
   fChainPhysics->SetBranchAddress("mu_nInnerSmallHits", &mu_nInnerSmallHits, &b_mu_nInnerSmallHits);
   fChainPhysics->SetBranchAddress("mu_nInnerLargeHits", &mu_nInnerLargeHits, &b_mu_nInnerLargeHits);
   fChainPhysics->SetBranchAddress("mu_nMiddleSmallHits", &mu_nMiddleSmallHits, &b_mu_nMiddleSmallHits);
   fChainPhysics->SetBranchAddress("mu_nMiddleLargeHits", &mu_nMiddleLargeHits, &b_mu_nMiddleLargeHits);
   fChainPhysics->SetBranchAddress("mu_nOuterSmallHits", &mu_nOuterSmallHits, &b_mu_nOuterSmallHits);
   fChainPhysics->SetBranchAddress("mu_nOuterLargeHits", &mu_nOuterLargeHits, &b_mu_nOuterLargeHits);
   fChainPhysics->SetBranchAddress("mu_nExtendedSmallHits", &mu_nExtendedSmallHits, &b_mu_nExtendedSmallHits);
   fChainPhysics->SetBranchAddress("mu_nExtendedLargeHits", &mu_nExtendedLargeHits, &b_mu_nExtendedLargeHits);
   fChainPhysics->SetBranchAddress("mu_nInnerSmallHoles", &mu_nInnerSmallHoles, &b_mu_nInnerSmallHoles);
   fChainPhysics->SetBranchAddress("mu_nInnerLargeHoles", &mu_nInnerLargeHoles, &b_mu_nInnerLargeHoles);
   fChainPhysics->SetBranchAddress("mu_nMiddleSmallHoles", &mu_nMiddleSmallHoles, &b_mu_nMiddleSmallHoles);
   fChainPhysics->SetBranchAddress("mu_nMiddleLargeHoles", &mu_nMiddleLargeHoles, &b_mu_nMiddleLargeHoles);
   fChainPhysics->SetBranchAddress("mu_nOuterSmallHoles", &mu_nOuterSmallHoles, &b_mu_nOuterSmallHoles);
   fChainPhysics->SetBranchAddress("mu_nOuterLargeHoles", &mu_nOuterLargeHoles, &b_mu_nOuterLargeHoles);
   fChainPhysics->SetBranchAddress("mu_nExtendedSmallHoles", &mu_nExtendedSmallHoles, &b_mu_nExtendedSmallHoles);
   fChainPhysics->SetBranchAddress("mu_nExtendedLargeHoles", &mu_nExtendedLargeHoles, &b_mu_nExtendedLargeHoles);
   fChainPhysics->SetBranchAddress("mu_nPhiLayer1Hits", &mu_nPhiLayer1Hits, &b_mu_nPhiLayer1Hits);
   fChainPhysics->SetBranchAddress("mu_nPhiLayer2Hits", &mu_nPhiLayer2Hits, &b_mu_nPhiLayer2Hits);
   fChainPhysics->SetBranchAddress("mu_nPhiLayer3Hits", &mu_nPhiLayer3Hits, &b_mu_nPhiLayer3Hits);
   fChainPhysics->SetBranchAddress("mu_nPhiLayer4Hits", &mu_nPhiLayer4Hits, &b_mu_nPhiLayer4Hits);
   fChainPhysics->SetBranchAddress("mu_nEtaLayer1Hits", &mu_nEtaLayer1Hits, &b_mu_nEtaLayer1Hits);
   fChainPhysics->SetBranchAddress("mu_nEtaLayer2Hits", &mu_nEtaLayer2Hits, &b_mu_nEtaLayer2Hits);
   fChainPhysics->SetBranchAddress("mu_nEtaLayer3Hits", &mu_nEtaLayer3Hits, &b_mu_nEtaLayer3Hits);
   fChainPhysics->SetBranchAddress("mu_nEtaLayer4Hits", &mu_nEtaLayer4Hits, &b_mu_nEtaLayer4Hits);
   fChainPhysics->SetBranchAddress("mu_nPhiLayer1Holes", &mu_nPhiLayer1Holes, &b_mu_nPhiLayer1Holes);
   fChainPhysics->SetBranchAddress("mu_nPhiLayer2Holes", &mu_nPhiLayer2Holes, &b_mu_nPhiLayer2Holes);
   fChainPhysics->SetBranchAddress("mu_nPhiLayer3Holes", &mu_nPhiLayer3Holes, &b_mu_nPhiLayer3Holes);
   fChainPhysics->SetBranchAddress("mu_nPhiLayer4Holes", &mu_nPhiLayer4Holes, &b_mu_nPhiLayer4Holes);
   fChainPhysics->SetBranchAddress("mu_nEtaLayer1Holes", &mu_nEtaLayer1Holes, &b_mu_nEtaLayer1Holes);
   fChainPhysics->SetBranchAddress("mu_nEtaLayer2Holes", &mu_nEtaLayer2Holes, &b_mu_nEtaLayer2Holes);
   fChainPhysics->SetBranchAddress("mu_nEtaLayer3Holes", &mu_nEtaLayer3Holes, &b_mu_nEtaLayer3Holes);
   fChainPhysics->SetBranchAddress("mu_nEtaLayer4Holes", &mu_nEtaLayer4Holes, &b_mu_nEtaLayer4Holes);
   fChainPhysics->SetBranchAddress("mu_cb_d0", &mu_cb_d0, &b_mu_cb_d0);
   fChainPhysics->SetBranchAddress("mu_cb_z0", &mu_cb_z0, &b_mu_cb_z0);
   fChainPhysics->SetBranchAddress("mu_cb_phi0", &mu_cb_phi0, &b_mu_cb_phi0);
   fChainPhysics->SetBranchAddress("mu_cb_theta", &mu_cb_theta, &b_mu_cb_theta);
   fChainPhysics->SetBranchAddress("mu_cb_qOverP", &mu_cb_qOverP, &b_mu_cb_qOverP);
   fChainPhysics->SetBranchAddress("mu_cb_vx", &mu_cb_vx, &b_mu_cb_vx);
   fChainPhysics->SetBranchAddress("mu_cb_vy", &mu_cb_vy, &b_mu_cb_vy);
   fChainPhysics->SetBranchAddress("mu_cb_vz", &mu_cb_vz, &b_mu_cb_vz);
   fChainPhysics->SetBranchAddress("museg_n", &museg_n, &b_TGC_museg_n);
   fChainPhysics->SetBranchAddress("museg_x", &museg_x, &b_museg_x);
   fChainPhysics->SetBranchAddress("museg_y", &museg_y, &b_museg_y);
   fChainPhysics->SetBranchAddress("museg_z", &museg_z, &b_museg_z);
   fChainPhysics->SetBranchAddress("museg_px", &museg_px, &b_museg_px);
   fChainPhysics->SetBranchAddress("museg_py", &museg_py, &b_museg_py);
   fChainPhysics->SetBranchAddress("museg_pz", &museg_pz, &b_museg_pz);
   fChainPhysics->SetBranchAddress("museg_t0", &museg_t0, &b_museg_t0);
   fChainPhysics->SetBranchAddress("museg_t0error", &museg_t0error, &b_museg_t0error);
   fChainPhysics->SetBranchAddress("museg_chi2", &museg_chi2, &b_museg_chi2);
   fChainPhysics->SetBranchAddress("museg_ndof", &museg_ndof, &b_museg_ndof);
   fChainPhysics->SetBranchAddress("museg_sector", &museg_sector, &b_museg_sector);
   fChainPhysics->SetBranchAddress("museg_stationName", &museg_stationName, &b_museg_stationName);
   fChainPhysics->SetBranchAddress("museg_stationEta", &museg_stationEta, &b_museg_stationEta);
   fChainPhysics->SetBranchAddress("museg_author", &museg_author, &b_museg_author);
   fChainPhysics->SetBranchAddress("ext_mu_bias_n", &ext_mu_bias_n, &b_ext_mu_bias_n);
   fChainPhysics->SetBranchAddress("ext_mu_bias_type", &ext_mu_bias_type, &b_ext_mu_bias_type);
   fChainPhysics->SetBranchAddress("ext_mu_bias_index", &ext_mu_bias_index, &b_ext_mu_bias_index);
   fChainPhysics->SetBranchAddress("ext_mu_bias_size", &ext_mu_bias_size, &b_ext_mu_bias_size);
   fChainPhysics->SetBranchAddress("ext_mu_bias_targetVec", &ext_mu_bias_targetVec, &b_ext_mu_bias_targetVec);
   fChainPhysics->SetBranchAddress("ext_mu_bias_targetDistanceVec", &ext_mu_bias_targetDistanceVec, &b_ext_mu_bias_targetDistanceVec);
   fChainPhysics->SetBranchAddress("ext_mu_bias_targetEtaVec", &ext_mu_bias_targetEtaVec, &b_ext_mu_bias_targetEtaVec);
   fChainPhysics->SetBranchAddress("ext_mu_bias_targetPhiVec", &ext_mu_bias_targetPhiVec, &b_ext_mu_bias_targetPhiVec);
   fChainPhysics->SetBranchAddress("ext_mu_bias_targetDeltaEtaVec", &ext_mu_bias_targetDeltaEtaVec, &b_ext_mu_bias_targetDeltaEtaVec);
   fChainPhysics->SetBranchAddress("ext_mu_bias_targetDeltaPhiVec", &ext_mu_bias_targetDeltaPhiVec, &b_ext_mu_bias_targetDeltaPhiVec);
   fChainPhysics->SetBranchAddress("ext_mu_bias_targetPxVec", &ext_mu_bias_targetPxVec, &b_ext_mu_bias_targetPxVec);
   fChainPhysics->SetBranchAddress("ext_mu_bias_targetPyVec", &ext_mu_bias_targetPyVec, &b_ext_mu_bias_targetPyVec);
   fChainPhysics->SetBranchAddress("ext_mu_bias_targetPzVec", &ext_mu_bias_targetPzVec, &b_ext_mu_bias_targetPzVec);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_n", &ext_mu_ubias_n, &b_ext_mu_ubias_n);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_type", &ext_mu_ubias_type, &b_ext_mu_ubias_type);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_index", &ext_mu_ubias_index, &b_ext_mu_ubias_index);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_size", &ext_mu_ubias_size, &b_ext_mu_ubias_size);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_targetVec", &ext_mu_ubias_targetVec, &b_ext_mu_ubias_targetVec);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_targetDistanceVec", &ext_mu_ubias_targetDistanceVec, &b_ext_mu_ubias_targetDistanceVec);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_targetEtaVec", &ext_mu_ubias_targetEtaVec, &b_ext_mu_ubias_targetEtaVec);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_targetPhiVec", &ext_mu_ubias_targetPhiVec, &b_ext_mu_ubias_targetPhiVec);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_targetDeltaEtaVec", &ext_mu_ubias_targetDeltaEtaVec, &b_ext_mu_ubias_targetDeltaEtaVec);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_targetDeltaPhiVec", &ext_mu_ubias_targetDeltaPhiVec, &b_ext_mu_ubias_targetDeltaPhiVec);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_targetPxVec", &ext_mu_ubias_targetPxVec, &b_ext_mu_ubias_targetPxVec);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_targetPyVec", &ext_mu_ubias_targetPyVec, &b_ext_mu_ubias_targetPyVec);
   fChainPhysics->SetBranchAddress("ext_mu_ubias_targetPzVec", &ext_mu_ubias_targetPzVec, &b_ext_mu_ubias_targetPzVec);
   fChainPhysics->SetBranchAddress("trigger_info_n", &trigger_info_n, &b_trigger_info_n);
   fChainPhysics->SetBranchAddress("trigger_info_chain", &trigger_info_chain, &b_trigger_info_chain);
   fChainPhysics->SetBranchAddress("trigger_info_isPassed", &trigger_info_isPassed, &b_trigger_info_isPassed);
   fChainPhysics->SetBranchAddress("trigger_info_nTracks", &trigger_info_nTracks, &b_trigger_info_nTracks);
   fChainPhysics->SetBranchAddress("trigger_info_typeVec", &trigger_info_typeVec, &b_trigger_info_typeVec);
   fChainPhysics->SetBranchAddress("trigger_info_ptVec", &trigger_info_ptVec, &b_trigger_info_ptVec);
   fChainPhysics->SetBranchAddress("trigger_info_etaVec", &trigger_info_etaVec, &b_trigger_info_etaVec);
   fChainPhysics->SetBranchAddress("trigger_info_phiVec", &trigger_info_phiVec, &b_trigger_info_phiVec);
   fChainPhysics->SetBranchAddress("vxp_n", &vxp_n, &b_vxp_n);
   fChainPhysics->SetBranchAddress("vxp_x", &vxp_x, &b_vxp_x);
   fChainPhysics->SetBranchAddress("vxp_y", &vxp_y, &b_vxp_y);
   fChainPhysics->SetBranchAddress("vxp_z", &vxp_z, &b_vxp_z);
   fChainPhysics->SetBranchAddress("vxp_cov_x", &vxp_cov_x, &b_vxp_cov_x);
   fChainPhysics->SetBranchAddress("vxp_cov_y", &vxp_cov_y, &b_vxp_cov_y);
   fChainPhysics->SetBranchAddress("vxp_cov_z", &vxp_cov_z, &b_vxp_cov_z);
   fChainPhysics->SetBranchAddress("vxp_cov_xy", &vxp_cov_xy, &b_vxp_cov_xy);
   fChainPhysics->SetBranchAddress("vxp_cov_xz", &vxp_cov_xz, &b_vxp_cov_xz);
   fChainPhysics->SetBranchAddress("vxp_cov_yz", &vxp_cov_yz, &b_vxp_cov_yz);
   fChainPhysics->SetBranchAddress("vxp_chi2", &vxp_chi2, &b_vxp_chi2);
   fChainPhysics->SetBranchAddress("vxp_ndof", &vxp_ndof, &b_vxp_ndof);
   fChainPhysics->SetBranchAddress("vxp_nTracks", &vxp_nTracks, &b_vxp_nTracks);
   fChainPhysics->SetBranchAddress("vxp_type", &vxp_type, &b_vxp_type);
   fChainPhysics->SetBranchAddress("TGC_prd_n", &TGC_prd_n, &b_TGC_prd_n);
   fChainPhysics->SetBranchAddress("TGC_prd_x", &TGC_prd_x, &b_TGC_prd_x);
   fChainPhysics->SetBranchAddress("TGC_prd_y", &TGC_prd_y, &b_TGC_prd_y);
   fChainPhysics->SetBranchAddress("TGC_prd_z", &TGC_prd_z, &b_TGC_prd_z);
   fChainPhysics->SetBranchAddress("TGC_prd_shortWidth", &TGC_prd_shortWidth, &b_TGC_prd_shortWidth);
   fChainPhysics->SetBranchAddress("TGC_prd_longWidth", &TGC_prd_longWidth, &b_TGC_prd_longWidth);
   fChainPhysics->SetBranchAddress("TGC_prd_length", &TGC_prd_length, &b_TGC_prd_length);
   fChainPhysics->SetBranchAddress("TGC_prd_isStrip", &TGC_prd_isStrip, &b_TGC_prd_isStrip);
   fChainPhysics->SetBranchAddress("TGC_prd_gasGap", &TGC_prd_gasGap, &b_TGC_prd_gasGap);
   fChainPhysics->SetBranchAddress("TGC_prd_channel", &TGC_prd_channel, &b_TGC_prd_channel);
   fChainPhysics->SetBranchAddress("TGC_prd_eta", &TGC_prd_eta, &b_TGC_prd_eta);
   fChainPhysics->SetBranchAddress("TGC_prd_phi", &TGC_prd_phi, &b_TGC_prd_phi);
   fChainPhysics->SetBranchAddress("TGC_prd_station", &TGC_prd_station, &b_TGC_prd_station);
   fChainPhysics->SetBranchAddress("TGC_prd_bunch", &TGC_prd_bunch, &b_TGC_prd_bunch);
   fChainPhysics->SetBranchAddress("MDT_prd_n", &MDT_prd_n, &b_MDT_prd_n);
   fChainPhysics->SetBranchAddress("MDT_prd_x", &MDT_prd_x, &b_MDT_prd_x);
   fChainPhysics->SetBranchAddress("MDT_prd_y", &MDT_prd_y, &b_MDT_prd_y);
   fChainPhysics->SetBranchAddress("MDT_prd_z", &MDT_prd_z, &b_MDT_prd_z);
   fChainPhysics->SetBranchAddress("MDT_prd_adc", &MDT_prd_adc, &b_MDT_prd_adc);
   fChainPhysics->SetBranchAddress("MDT_prd_tdc", &MDT_prd_tdc, &b_MDT_prd_tdc);
   fChainPhysics->SetBranchAddress("MDT_prd_status", &MDT_prd_status, &b_MDT_prd_status);
   fChainPhysics->SetBranchAddress("MDT_prd_drift_radius", &MDT_prd_drift_radius, &b_MDT_prd_drift_radius);
   fChainPhysics->SetBranchAddress("MDT_prd_drift_radius_error", &MDT_prd_drift_radius_error, &b_MDT_prd_drift_radius_error);
   fChainPhysics->SetBranchAddress("TGC_coin_n", &TGC_coin_n, &b_TGC_coin_n);
   fChainPhysics->SetBranchAddress("TGC_coin_x_In", &TGC_coin_x_In, &b_TGC_coin_x_In);
   fChainPhysics->SetBranchAddress("TGC_coin_y_In", &TGC_coin_y_In, &b_TGC_coin_y_In);
   fChainPhysics->SetBranchAddress("TGC_coin_z_In", &TGC_coin_z_In, &b_TGC_coin_z_In);
   fChainPhysics->SetBranchAddress("TGC_coin_x_Out", &TGC_coin_x_Out, &b_TGC_coin_x_Out);
   fChainPhysics->SetBranchAddress("TGC_coin_y_Out", &TGC_coin_y_Out, &b_TGC_coin_y_Out);
   fChainPhysics->SetBranchAddress("TGC_coin_z_Out", &TGC_coin_z_Out, &b_TGC_coin_z_Out);
   fChainPhysics->SetBranchAddress("TGC_coin_width_In", &TGC_coin_width_In, &b_TGC_coin_width_In);
   fChainPhysics->SetBranchAddress("TGC_coin_width_Out", &TGC_coin_width_Out, &b_TGC_coin_width_Out);
   fChainPhysics->SetBranchAddress("TGC_coin_width_R", &TGC_coin_width_R, &b_TGC_coin_width_R);
   fChainPhysics->SetBranchAddress("TGC_coin_width_Phi", &TGC_coin_width_Phi, &b_TGC_coin_width_Phi);
   fChainPhysics->SetBranchAddress("TGC_coin_isAside", &TGC_coin_isAside, &b_TGC_coin_isAside);
   fChainPhysics->SetBranchAddress("TGC_coin_isForward", &TGC_coin_isForward, &b_TGC_coin_isForward);
   fChainPhysics->SetBranchAddress("TGC_coin_isStrip", &TGC_coin_isStrip, &b_TGC_coin_isStrip);
   fChainPhysics->SetBranchAddress("TGC_coin_isInner", &TGC_coin_isInner, &b_TGC_coin_isInner);
   fChainPhysics->SetBranchAddress("TGC_coin_isPositiveDeltaR", &TGC_coin_isPositiveDeltaR, &b_TGC_coin_isPositiveDeltaR);
   fChainPhysics->SetBranchAddress("TGC_coin_type", &TGC_coin_type, &b_TGC_coin_type);
   fChainPhysics->SetBranchAddress("TGC_coin_trackletId", &TGC_coin_trackletId, &b_TGC_coin_trackletId);
   fChainPhysics->SetBranchAddress("TGC_coin_trackletIdStrip", &TGC_coin_trackletIdStrip, &b_TGC_coin_trackletIdStrip);
   fChainPhysics->SetBranchAddress("TGC_coin_phi", &TGC_coin_phi, &b_TGC_coin_phi);
   fChainPhysics->SetBranchAddress("TGC_coin_roi", &TGC_coin_roi, &b_TGC_coin_roi);
   fChainPhysics->SetBranchAddress("TGC_coin_pt", &TGC_coin_pt, &b_TGC_coin_pt);
   fChainPhysics->SetBranchAddress("TGC_coin_delta", &TGC_coin_delta, &b_TGC_coin_delta);
   fChainPhysics->SetBranchAddress("TGC_coin_sub", &TGC_coin_sub, &b_TGC_coin_sub);
   fChainPhysics->SetBranchAddress("TGC_coin_veto", &TGC_coin_veto, &b_TGC_coin_veto);
   fChainPhysics->SetBranchAddress("TGC_coin_bunch", &TGC_coin_bunch, &b_TGC_coin_bunch);
   fChainPhysics->SetBranchAddress("TGC_coin_inner", &TGC_coin_inner, &b_TGC_coin_inner);
   fChainPhysics->SetBranchAddress("TILE_murcv_trig_n", &TILE_murcv_trig_n, &b_TILE_murcv_trig_n);
   fChainPhysics->SetBranchAddress("TILE_murcv_trig_mod", &TILE_murcv_trig_mod, &b_TILE_murcv_trig_mod);
   fChainPhysics->SetBranchAddress("TILE_murcv_trig_part", &TILE_murcv_trig_part, &b_TILE_murcv_trig_part);
   fChainPhysics->SetBranchAddress("TILE_murcv_trig_bit0", &TILE_murcv_trig_bit0, &b_TILE_murcv_trig_bit0);
   fChainPhysics->SetBranchAddress("TILE_murcv_trig_bit1", &TILE_murcv_trig_bit1, &b_TILE_murcv_trig_bit1);
   fChainPhysics->SetBranchAddress("TILE_murcv_trig_bit2", &TILE_murcv_trig_bit2, &b_TILE_murcv_trig_bit2);
   fChainPhysics->SetBranchAddress("TILE_murcv_trig_bit3", &TILE_murcv_trig_bit3, &b_TILE_murcv_trig_bit3);
   fChainPhysics->SetBranchAddress("TILE_murcv_raw_n", &TILE_murcv_raw_n, &b_TILE_murcv_raw_n);
   fChainPhysics->SetBranchAddress("TILE_murcv_raw_count", &TILE_murcv_raw_count, &b_TILE_murcv_raw_count);
   fChainPhysics->SetBranchAddress("TILE_murcv_raw_energy", &TILE_murcv_raw_energy, &b_TILE_murcv_raw_energy);
   fChainPhysics->SetBranchAddress("TILE_murcv_raw_ros", &TILE_murcv_raw_ros, &b_TILE_murcv_raw_ros);
   fChainPhysics->SetBranchAddress("TILE_murcv_raw_drawer", &TILE_murcv_raw_drawer, &b_TILE_murcv_raw_drawer);
   fChainPhysics->SetBranchAddress("TILE_murcv_raw_channel", &TILE_murcv_raw_channel, &b_TILE_murcv_raw_channel);
   fChainPhysics->SetBranchAddress("TILE_murcv_digit_n", &TILE_murcv_digit_n, &b_TILE_murcv_digit_n);
   fChainPhysics->SetBranchAddress("TILE_murcv_digit_nSamples", &TILE_murcv_digit_nSamples, &b_TILE_murcv_digit_nSamples);
   fChainPhysics->SetBranchAddress("TILE_murcv_digit_ros", &TILE_murcv_digit_ros, &b_TILE_murcv_digit_ros);
   fChainPhysics->SetBranchAddress("TILE_murcv_digit_drawer", &TILE_murcv_digit_drawer, &b_TILE_murcv_digit_drawer);
   fChainPhysics->SetBranchAddress("TILE_murcv_digit_channel", &TILE_murcv_digit_channel, &b_TILE_murcv_digit_channel);
   fChainPhysics->SetBranchAddress("TILE_murcv_digit_sampleVec", &TILE_murcv_digit_sampleVec, &b_TILE_murcv_digit_sampleVec);
   fChainPhysics->SetBranchAddress("TGC_hierarchy_n", &TGC_hierarchy_n, &b_TGC_hierarchy_n);
   fChainPhysics->SetBranchAddress("TGC_hierarchy_index", &TGC_hierarchy_index, &b_TGC_hierarchy_index);
   fChainPhysics->SetBranchAddress("TGC_hierarchy_dR_hiPt", &TGC_hierarchy_dR_hiPt, &b_TGC_hierarchy_dR_hiPt);
   fChainPhysics->SetBranchAddress("TGC_hierarchy_dPhi_hiPt", &TGC_hierarchy_dPhi_hiPt, &b_TGC_hierarchy_dPhi_hiPt);
   fChainPhysics->SetBranchAddress("TGC_hierarchy_dR_tracklet", &TGC_hierarchy_dR_tracklet, &b_TGC_hierarchy_dR_tracklet);
   fChainPhysics->SetBranchAddress("TGC_hierarchy_dPhi_tracklet", &TGC_hierarchy_dPhi_tracklet, &b_TGC_hierarchy_dPhi_tracklet);
   fChainPhysics->SetBranchAddress("TGC_hierarchy_isChamberBoundary", &TGC_hierarchy_isChamberBoundary, &b_TGC_hierarchy_isChamberBoundary);
   fChainPhysics->SetBranchAddress("muctpi_candidateMultiplicities", &muctpi_candidateMultiplicities, &b_muctpi_candidateMultiplicities);
   fChainPhysics->SetBranchAddress("muctpi_nDataWords", &muctpi_nDataWords, &b_muctpi_nDataWords);
   fChainPhysics->SetBranchAddress("muctpi_dataWords", &muctpi_dataWords, &b_muctpi_dataWords);
   fChainPhysics->SetBranchAddress("muctpi_dw_eta", &muctpi_dw_eta, &b_muctpi_dw_eta);
   fChainPhysics->SetBranchAddress("muctpi_dw_phi", &muctpi_dw_phi, &b_muctpi_dw_phi);
   fChainPhysics->SetBranchAddress("muctpi_dw_source", &muctpi_dw_source, &b_muctpi_dw_source);
   fChainPhysics->SetBranchAddress("muctpi_dw_hemisphere", &muctpi_dw_hemisphere, &b_muctpi_dw_hemisphere);
   fChainPhysics->SetBranchAddress("muctpi_dw_bcid", &muctpi_dw_bcid, &b_muctpi_dw_bcid);
   fChainPhysics->SetBranchAddress("muctpi_dw_sectorID", &muctpi_dw_sectorID, &b_muctpi_dw_sectorID);
   fChainPhysics->SetBranchAddress("muctpi_dw_thrNumber", &muctpi_dw_thrNumber, &b_muctpi_dw_thrNumber);
   fChainPhysics->SetBranchAddress("muctpi_dw_roi", &muctpi_dw_roi, &b_muctpi_dw_roi);
   fChainPhysics->SetBranchAddress("muctpi_dw_veto", &muctpi_dw_veto, &b_muctpi_dw_veto);
   fChainPhysics->SetBranchAddress("muctpi_dw_firstCandidate", &muctpi_dw_firstCandidate, &b_muctpi_dw_firstCandidate);
   fChainPhysics->SetBranchAddress("muctpi_dw_moreCandInRoI", &muctpi_dw_moreCandInRoI, &b_muctpi_dw_moreCandInRoI);
   fChainPhysics->SetBranchAddress("muctpi_dw_moreCandInSector", &muctpi_dw_moreCandInSector, &b_muctpi_dw_moreCandInSector);
   fChainPhysics->SetBranchAddress("muctpi_dw_charge", &muctpi_dw_charge, &b_muctpi_dw_charge);
   fChainPhysics->SetBranchAddress("muctpi_dw_candidateVetoed", &muctpi_dw_candidateVetoed, &b_muctpi_dw_candidateVetoed);
  //  fChainPhysics->SetBranchAddress("TILE_cell_n", &TILE_cell_n, &b_TILE_cell_n);
  //  fChainPhysics->SetBranchAddress("TILE_cell_E", &TILE_cell_E, &b_TILE_cell_E);
  //  fChainPhysics->SetBranchAddress("TILE_cell_Et", &TILE_cell_Et, &b_TILE_cell_Et);
  //  fChainPhysics->SetBranchAddress("TILE_cell_eta", &TILE_cell_eta, &b_TILE_cell_eta);
  //  fChainPhysics->SetBranchAddress("TILE_cell_phi", &TILE_cell_phi, &b_TILE_cell_phi);
  //  fChainPhysics->SetBranchAddress("TILE_cell_sinTh", &TILE_cell_sinTh, &b_TILE_cell_sinTh);
  //  fChainPhysics->SetBranchAddress("TILE_cell_cosTh", &TILE_cell_cosTh, &b_TILE_cell_cosTh);
  //  fChainPhysics->SetBranchAddress("TILE_cell_cotTh", &TILE_cell_cotTh, &b_TILE_cell_cotTh);
  //  fChainPhysics->SetBranchAddress("TILE_cell_x", &TILE_cell_x, &b_TILE_cell_x);
  //  fChainPhysics->SetBranchAddress("TILE_cell_y", &TILE_cell_y, &b_TILE_cell_y);
  //  fChainPhysics->SetBranchAddress("TILE_cell_z", &TILE_cell_z, &b_TILE_cell_z);
  //  fChainPhysics->SetBranchAddress("TILE_cell_badcell", &TILE_cell_badcell, &b_TILE_cell_badcell);
  //  fChainPhysics->SetBranchAddress("TILE_cell_partition", &TILE_cell_partition, &b_TILE_cell_partition);
  //  fChainPhysics->SetBranchAddress("TILE_cell_section", &TILE_cell_section, &b_TILE_cell_section);
  //  fChainPhysics->SetBranchAddress("TILE_cell_side", &TILE_cell_side, &b_TILE_cell_side);
  //  fChainPhysics->SetBranchAddress("TILE_cell_module", &TILE_cell_module, &b_TILE_cell_module);
  //  fChainPhysics->SetBranchAddress("TILE_cell_tower", &TILE_cell_tower, &b_TILE_cell_tower);
  //  fChainPhysics->SetBranchAddress("TILE_cell_sample", &TILE_cell_sample, &b_TILE_cell_sample);

   Notify();
}

Bool_t TMDB::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TMDB::ShowPhysics(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChainPhysics) return;
   fChainPhysics->Show(entry);
}

Int_t TMDB::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TMDB_cxx
