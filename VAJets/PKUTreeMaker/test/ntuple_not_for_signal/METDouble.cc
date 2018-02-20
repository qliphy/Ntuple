// -*- C++ -*-
//
// Package:    METDouble
// Class:      METDouble
// 
/**\class METDouble METDouble.cc RA2Classic/METDouble/src/METDouble.cc
 * 
 * Description: [one line class summary]
 * 
 * Implementation:
 *     [Notes on implementation]
 */
//
// Original Author:  Arne-Rasmus Draeger,68/111,4719,
//         Created:  Fri Apr 11 16:35:33 CEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrectorCalculator.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
//
// class declaration
//

class METDouble : public edm::EDProducer {
	public:
		explicit METDouble(const edm::ParameterSet&);
		~METDouble();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() ;
		virtual void produce(edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run&, edm::EventSetup const&);
		virtual void endRun(edm::Run&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

		edm::EDGetTokenT<edm::View<pat::MET> > metToken_;  
		edm::EDGetTokenT<edm::View<pat::Jet> > JetToken_;
		edm::EDGetTokenT<edm::View<pat::Muon> > MuToken_;
		edm::EDGetTokenT<double> RhoToken_;

		//  std::vector<std::string> jecPayloadNames_;
		bool corrMet;
		double corrEx;
		double corrEy;
		bool corrMet_up;
		double corrEx_up;
		double corrEy_up;
		bool corrMet_down;
		double corrEx_down;
		double corrEy_down;
		double corrSumEt;
		std::string l1file;
		std::string l2file;
		std::string l3file;
		std::string l2l3file;
		std::string uncfile;
		bool doJEC;
		JetCorrectorParameters *L2L3JetPar;

		// ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
METDouble::METDouble(const edm::ParameterSet& iConfig):
	metToken_(consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag> ("METTag"))),
	JetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("JetTag"))),
	MuToken_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag> ("MuTag"))),
	RhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("RhoTag")))
	//  jecPayloadNames_( iConfig.getParameter<std::vector<std::string> >("jecPayloadNames") ) // JEC level payloads
{
	//register your produc
	corrMet = iConfig.getParameter<bool> ("corrMet");
	doJEC = iConfig.getParameter<bool> ("doJEC");
	l1file = iConfig.getParameter<std::string> ("L1File");
	l2file = iConfig.getParameter<std::string> ("L2File");
	l3file = iConfig.getParameter<std::string> ("L3File");
	l2l3file = iConfig.getParameter<std::string> ("L2L3File");
	uncfile = iConfig.getParameter<std::string>  ("uncFile");

	produces<double>("Pt");
	produces<double>("Phi");
	produces<double>("PtRaw");
	produces<double>("PhiRaw");
	produces<double>("CaloMetPt");
	produces<double>("CaloMetPhi");
	produces<double>("PtUp");
	produces<double>("PhiUp");
	produces<double>("PtDown");
	produces<double>("PhiDown");
	/* Examples
	 *   produces<ExampleData2>();
	 * 
	 *   //if do put with a label
	 *   produces<ExampleData2>("label");
	 * 
	 *   //if you want to put into the Run
	 *   produces<ExampleData2,InRun>();
	 */
	//now do what ever other initialization is needed

}


METDouble::~METDouble()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
	void
METDouble::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	double metpt_=0, metphi_=0;
	double metpt_up_=0, metphi_up_=0;
	double metpt_down_=0, metphi_down_=0;
	double rawmetpt_=0, rawmetphi_=0;
	double calometpt_=0, calometphi_=0;
	edm::Handle< edm::View<pat::MET> > MET;
	iEvent.getByToken(metToken_,MET); 
	edm::Handle<double> rho_ ;
	iEvent.getByToken(RhoToken_, rho_);
	edm::Handle< edm::View<pat::Jet> > Jets;
	iEvent.getByToken(JetToken_,Jets);
	edm::Handle< edm::View<pat::Muon> > Muons;
	iEvent.getByToken(MuToken_,Muons);

	//  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
	std::vector<JetCorrectorParameters> vPar;
	/*	for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNames_.begin(),
		payloadEnd = jecPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
		JetCorrectorParameters pars(*ipayload);
		vPar.push_back(pars);
		}*/

	std::vector<JetCorrectorParameters> vParL1;
	//	vParL1.push_back(JetCorrectorParameters(jecPayloadNames_[0]));

	if (l2l3file!="NONE")
		L2L3JetPar  = new JetCorrectorParameters(l2l3file);
	JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(l3file);
	JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(l2file);
	JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(l1file);

	vPar.push_back(*L1JetPar);
	vPar.push_back(*L2JetPar);
	vPar.push_back(*L3JetPar);
	if (l2l3file!="NONE")
		vPar.push_back(*L2L3JetPar);
	vParL1.push_back(*L1JetPar);

	FactorizedJetCorrector *JetCorrector = new FactorizedJetCorrector(vPar);
	FactorizedJetCorrector *JetCorrectorL1 = new FactorizedJetCorrector(vParL1);

	JetCorrectionUncertainty *JetUnc  = new JetCorrectionUncertainty(uncfile);

	double corrUp;
	double corrDown;

	if (!doJEC)  corrMet=false; //it does not have any sense to correct met if no JEC applied

	if (corrMet) {
		corrEx = 0;
		corrEy = 0;
		corrSumEt = 0;
		corrEx_up = 0;
		corrEy_up = 0;
		corrEx_down = 0;
		corrEy_down = 0;

		bool skipEM_ = true;
		double skipEMfractionThreshold_ = 0.9;
		bool skipMuons_ = true;
		double jetCorrEtaMax_ = 9.9;
		double type1JetPtThreshold_ = 15.0;
		std::string skipMuonSelection_string = "isGlobalMuon | isStandAloneMuon";
		StringCutObjectSelector<reco::Candidate>* skipMuonSelection_ = new StringCutObjectSelector<reco::Candidate>(skipMuonSelection_string,true);

		edm::View<pat::Jet>::const_iterator ijet = Jets->begin()-1;
		for (const pat::Jet &jet : *Jets) {

			ijet++;

			reco::Candidate::LorentzVector uncorrJet;
			// The pat::Jet "knows" if it has been corrected, so here
			// we can "uncorrect" the entire jet to apply the corrections
			// we want here.
			pat::Jet const * pJet = dynamic_cast<pat::Jet const *>( &*ijet );
			if ( pJet != 0 ) {
				uncorrJet = pJet->correctedP4(0);
			}
			// Otherwise, if we do not have pat::Jets on input, we just assume
			// the user has not corrected them upstream and use it as raw.
			else {
				uncorrJet = ijet->p4();
			}

			if (fabs(uncorrJet.eta()) < jetCorrEtaMax_)
				JetCorrector->setJetEta(uncorrJet.eta());
			else
				JetCorrector->setJetEta(TMath::Sign(1.,uncorrJet.eta())*jetCorrEtaMax_);
			JetCorrector->setJetPt(uncorrJet.pt());
			JetCorrector->setJetA(ijet->jetArea());
			JetCorrector->setRho(*(rho_.product())); 

			double corr = JetCorrector->getCorrection();

			JetUnc->setJetEta( uncorrJet.eta() );
			JetUnc->setJetPt( corr * uncorrJet.pt() );
			corrUp = corr * (1 + fabs(JetUnc->getUncertainty(1)));
			JetUnc->setJetEta( uncorrJet.eta() );
			JetUnc->setJetPt( corr * uncorrJet.pt() );
			corrDown = corr * ( 1 - fabs(JetUnc->getUncertainty(-1)) );

			double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
			if ( skipEM_ && emEnergyFraction > skipEMfractionThreshold_ ) continue;

			reco::Candidate::LorentzVector rawJetP4 = jet.correctedP4(0);

			//	    if ( skipMuons_ && jet.muonMultiplicity() != 0 ) {
			if ( skipMuons_ ) {
				const std::vector<reco::CandidatePtr> & cands = jet.daughterPtrVector();
				for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();
						cand != cands.end(); ++cand ) {
					const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
					const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
					if ( mu != 0 &&  (*skipMuonSelection_)(*mu) ) {
						reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
						rawJetP4 -= muonP4;
					}
				}
			}

			reco::Candidate::LorentzVector corrJetP4 = corr*rawJetP4;
			reco::Candidate::LorentzVector corrJetP4_up = corrUp*rawJetP4;
			reco::Candidate::LorentzVector corrJetP4_down = corrDown*rawJetP4;

			if ( corrJetP4.pt() > type1JetPtThreshold_) {
				reco::Candidate::LorentzVector tmpP4 = jet.correctedP4(0);
				if (fabs(tmpP4.eta()) < jetCorrEtaMax_)
					JetCorrectorL1->setJetEta(tmpP4.eta());
				else
					JetCorrectorL1->setJetEta(TMath::Sign(1.,tmpP4.eta())*jetCorrEtaMax_);
				JetCorrectorL1->setJetPt(tmpP4.pt());
				JetCorrectorL1->setJetA(ijet->jetArea());
				JetCorrectorL1->setRho(*(rho_.product())); 

				corr = JetCorrectorL1->getCorrection();

				if (fabs(tmpP4.eta()) < jetCorrEtaMax_)
					JetUnc->setJetEta( tmpP4.eta() );
				else
					JetUnc->setJetEta(TMath::Sign(1.,tmpP4.eta())*jetCorrEtaMax_);
				JetUnc->setJetPt( corr * tmpP4.pt() );
				corrUp = corr * (1 + fabs(JetUnc->getUncertainty(1)));

				if (fabs(tmpP4.eta()) < jetCorrEtaMax_)
					JetUnc->setJetEta( tmpP4.eta() );
				else
					JetUnc->setJetEta(TMath::Sign(1.,tmpP4.eta())*jetCorrEtaMax_);
				JetUnc->setJetPt( corr * tmpP4.pt() );
				corrDown = corr * (1 - fabs(JetUnc->getUncertainty(-1)));

				reco::Candidate::LorentzVector rawJetP4offsetCorr = corr*rawJetP4;
				corrEx -= (corrJetP4.px() - rawJetP4offsetCorr.px());
				corrEy -= (corrJetP4.py() - rawJetP4offsetCorr.py());
				corrSumEt += (corrJetP4.Et() - rawJetP4offsetCorr.Et());

				reco::Candidate::LorentzVector rawJetP4offsetCorr_up = corrUp*rawJetP4;
				corrEx_up -= (corrJetP4_up.px() - rawJetP4offsetCorr_up.px());
				corrEy_up -= (corrJetP4_up.py() - rawJetP4offsetCorr_up.py());

				reco::Candidate::LorentzVector rawJetP4offsetCorr_down = corrDown*rawJetP4;
				corrEx_down -= (corrJetP4_down.px() - rawJetP4offsetCorr_down.px());
				corrEy_down -= (corrJetP4_down.py() - rawJetP4offsetCorr_down.py());
			}

		}
		//	  pat::MET pmet = pat::MET( &*(MET->at(0)) );

		//const float rawPt = MET->at(0).shiftedPt(pat::MET::NoShift, pat::MET::Raw);
		//const float rawPhi = MET->at(0).shiftedPhi(pat::MET::NoShift, pat::MET::Raw);
		//	  const float rawSumEt = MET->at(0).shiftedSumEt(pat::MET::NoShift, pat::MET::Raw);
		const float rawPt = MET->at(0).uncorPt();//met.shiftedPt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
		const float rawPhi   = MET->at(0).uncorPhi();//met.shiftedPhi(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
		//	  const float rawSumEt = met.uncorSumEt();//met.shiftedSumEt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);

		TVector2 rawMET_;
		rawMET_.SetMagPhi (rawPt, rawPhi );
		Double_t rawPx = rawMET_.Px();
		Double_t rawPy = rawMET_.Py();
		//	  Double_t rawEt = std::hypot(rawPx,rawPy);

		double pxcorr = rawPx+corrEx;
		double pycorr = rawPy+corrEy;
		double et = std::hypot(pxcorr,pycorr);
		//	  double sumEtcorr = rawSumEt+corrSumEt;
		TLorentzVector corrmet; corrmet.SetPxPyPzE(pxcorr,pycorr,0.,et);

		metpt_= et;
		metphi_ = corrmet.Phi();
		rawmetpt_ = rawPt;
		rawmetphi_ = rawPhi;


		double pxcorr_up = rawPx+corrEx_up;
		double pycorr_up = rawPy+corrEy_up;
		double et_up = std::hypot(pxcorr_up,pycorr_up);
		//	  double sumEtcorr = rawSumEt+corrSumEt;
		TLorentzVector corrmet_up; corrmet_up.SetPxPyPzE(pxcorr_up,pycorr_up,0.,et_up);

		metpt_up_= et_up;
		metphi_up_ = corrmet_up.Phi();


		double pxcorr_down = rawPx+corrEx_down;
		double pycorr_down = rawPy+corrEy_down;
		double et_down = std::hypot(pxcorr_down,pycorr_down);
		//	  double sumEtcorr = rawSumEt+corrSumEt;
		TLorentzVector corrmet_down; corrmet_down.SetPxPyPzE(pxcorr_down,pycorr_down,0.,et_down);

		metpt_down_= et_down;
		metphi_down_ = corrmet_down.Phi();

		//	  std::cout<<"corrEx: "<<corrEx<<" rawPx: "<<rawPx<<" raw met: "<<rawPt<<" new met: "<<et<<std::endl;
		//metcorrPx_ = corrEx;
		//metCorrPy_ = corrEy;
		//metSumEt_ = sumEtcorr;
		}

		else {
			if(MET.isValid() ){
				metpt_=MET->at(0).pt();
				metphi_=MET->at(0).phi();
				rawmetpt_ = MET->at(0).pt();
				rawmetphi_ = MET->at(0).phi();
				calometpt_ = MET->at(0).caloMETPt();
				calometphi_ = MET->at(0).caloMETPhi();
			}	
			else std::cout<<"METDouble::Invlide Tag: "<<std::endl;
		}

		//	std::cout<<"met: "<<metpt_<<std::endl;

		delete JetCorrector;
		delete JetCorrectorL1;
		delete JetUnc;
		delete L1JetPar;
		delete L2JetPar;
		delete L3JetPar;

		std::auto_ptr<double> htp(new double(metpt_));
		iEvent.put(htp,"Pt");
		std::auto_ptr<double> htp2(new double(metphi_));
		iEvent.put(htp2,"Phi");
		std::auto_ptr<double> htp3(new double(rawmetpt_));
		iEvent.put(htp3,"PtRaw");
		std::auto_ptr<double> htp4(new double(rawmetphi_));
		iEvent.put(htp4,"PhiRaw");
		std::auto_ptr<double> htp5(new double(calometpt_));
		iEvent.put(htp5,"CaloMetPt");
		std::auto_ptr<double> htp6(new double(calometphi_));
		iEvent.put(htp6,"CaloMetPhi");
		std::auto_ptr<double> htp7(new double(metpt_up_));
		iEvent.put(htp7,"PtUp");
		std::auto_ptr<double> htp8(new double(metphi_up_));
		iEvent.put(htp8,"PhiUp");
		std::auto_ptr<double> htp9(new double(metpt_down_));
		iEvent.put(htp9,"PtDown");
		std::auto_ptr<double> htp10(new double(metphi_down_));
		iEvent.put(htp10,"PhiDown");

	}

	// ------------ method called once each job just before starting event loop  ------------
	void 
		METDouble::beginJob()
		{
		}

	// ------------ method called once each job just after ending the event loop  ------------
	void 
		METDouble::endJob() {
		}

	// ------------ method called when starting to processes a run  ------------
	void 
		METDouble::beginRun(edm::Run&, edm::EventSetup const&)
		{
		}

	// ------------ method called when ending the processing of a run  ------------
	void 
		METDouble::endRun(edm::Run&, edm::EventSetup const&)
		{
		}

	// ------------ method called when starting to processes a luminosity block  ------------
	void 
		METDouble::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
		{
		}

	// ------------ method called when ending the processing of a luminosity block  ------------
	void 
		METDouble::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
		{
		}

	// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
	void
		METDouble::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
			//The following says we do not know what parameters are allowed so do no validation
			// Please change this to state exactly what you do use, even if it is no parameters
			edm::ParameterSetDescription desc;
			desc.setUnknown();
			descriptions.addDefault(desc);
		}

	//define this as a plug-in
	DEFINE_FWK_MODULE(METDouble);
