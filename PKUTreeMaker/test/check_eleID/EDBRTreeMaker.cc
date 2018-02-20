// system include files
#include <iostream>
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"  
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include <TFormula.h>

#define Pi 3.141593
using namespace std;
class EDBRTreeMaker : public edm::EDAnalyzer {
public:
  explicit EDBRTreeMaker(const edm::ParameterSet&);
  ~EDBRTreeMaker();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
  virtual void endRun(const edm::Run&, const edm::EventSetup&) override;


  // ----------member data ---------------------------
  TTree* outTree_;
  double gen_ele_pt, gen_ele_eta, gen_ele_phi;
  double gen_tight_match_ele_pt, gen_tight_match_ele_eta, gen_tight_match_ele_phi, gen_tight_match_ele_dr;
  double gen_medium_match_ele_pt, gen_medium_match_ele_eta, gen_medium_match_ele_phi, gen_medium_match_ele_dr;
  double reco_ele_pt, reco_ele_eta, reco_ele_phi, reco_ele_dr;
  double reco_tight_match_ele_pt, reco_tight_match_ele_eta, reco_tight_match_ele_phi, reco_tight_match_ele_dr;
  double reco_medium_match_ele_pt, reco_medium_match_ele_eta, reco_medium_match_ele_phi, reco_medium_match_ele_dr;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genSrc_;
  edm::EDGetTokenT<edm::View<pat::Electron> > EleSrc2_;
  edm::EDGetTokenT<edm::View<pat::Electron> > EleSrc1_;
  edm::EDGetTokenT<edm::View<pat::Electron> > EleSrc0_;

  void setDummyValues();

};

//
// constructors and destructor
//
EDBRTreeMaker::EDBRTreeMaker(const edm::ParameterSet& iConfig)
{
  EleSrc2_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("EleSrc2")));
  EleSrc1_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("EleSrc1")));
  EleSrc0_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("EleSrc0")));
  genSrc_      = consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>( "genSrc") ) ;
  
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  outTree_ = fs->make<TTree>("EDBRCandidates","EDBR Candidates");
  /// Basic event quantities
  outTree_->Branch("gen_ele_pt"        ,&gen_ele_pt       ,"gen_ele_pt/D"       );
  outTree_->Branch("gen_ele_eta"        ,&gen_ele_eta       ,"gen_ele_eta/D"       );
  outTree_->Branch("gen_ele_phi"        ,&gen_ele_phi       ,"gen_ele_phi/D"       );
  outTree_->Branch("gen_tight_match_ele_pt"        ,&gen_tight_match_ele_pt       ,"gen_tight_match_ele_pt/D"       );
  outTree_->Branch("gen_tight_match_ele_eta"       ,&gen_tight_match_ele_eta      ,"gen_tight_match_ele_eta/D"       );
  outTree_->Branch("gen_tight_match_ele_phi"       ,&gen_tight_match_ele_phi     ,"gen_tight_match_ele_phi/D"       );
  outTree_->Branch("gen_tight_match_ele_dr"        ,&gen_tight_match_ele_dr       ,"gen_tight_match_ele_dr/D"       );
  outTree_->Branch("gen_medium_match_ele_pt"        ,&gen_medium_match_ele_pt       ,"gen_medium_match_ele_pt/D"       );
  outTree_->Branch("gen_medium_match_ele_eta"       ,&gen_medium_match_ele_eta      ,"gen_medium_match_ele_eta/D"       );
  outTree_->Branch("gen_medium_match_ele_phi"       ,&gen_medium_match_ele_phi     ,"gen_medium_match_ele_phi/D"       );
  outTree_->Branch("gen_medium_match_ele_dr"        ,&gen_medium_match_ele_dr       ,"gen_medium_match_ele_dr/D"       );
  outTree_->Branch("reco_ele_pt"        ,&reco_ele_pt       ,"reco_ele_pt/D"       );
  outTree_->Branch("reco_ele_eta"        ,&reco_ele_eta       ,"reco_ele_eta/D"       );
  outTree_->Branch("reco_ele_phi"        ,&reco_ele_phi       ,"reco_ele_phi/D"       );
  outTree_->Branch("reco_ele_dr"        ,&reco_ele_dr       ,"reco_ele_dr/D"       );
  outTree_->Branch("reco_tight_match_ele_pt"        ,&reco_tight_match_ele_pt       ,"reco_tight_match_ele_pt/D"       );
  outTree_->Branch("reco_tight_match_ele_eta"        ,&reco_tight_match_ele_eta       ,"reco_tight_match_ele_eta/D"       );
  outTree_->Branch("reco_tight_match_ele_phi"        ,&reco_tight_match_ele_phi       ,"reco_tight_match_ele_phi/D"       );
  outTree_->Branch("reco_tight_match_ele_dr"        ,&reco_tight_match_ele_dr       ,"reco_tight_match_ele_dr/D"       );
  outTree_->Branch("reco_medium_match_ele_pt"        ,&reco_medium_match_ele_pt       ,"reco_medium_match_ele_pt/D"       );
  outTree_->Branch("reco_medium_match_ele_eta"        ,&reco_medium_match_ele_eta       ,"reco_medium_match_ele_eta/D"       );
  outTree_->Branch("reco_medium_match_ele_phi"        ,&reco_medium_match_ele_phi       ,"reco_medium_match_ele_phi/D"       );
  outTree_->Branch("reco_medium_match_ele_dr"        ,&reco_medium_match_ele_dr       ,"reco_medium_match_ele_dr/D"       );
}


EDBRTreeMaker::~EDBRTreeMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


// ------------ method called for each event  ------------
void
EDBRTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   setDummyValues(); //Initalize variables with dummy values


   edm::Handle<edm::View<reco::GenParticle> > genParticles;//define genParticle
   iEvent.getByToken(genSrc_, genParticles);

   edm::Handle<edm::View<pat::Electron>> eles_tight;
   iEvent.getByToken(EleSrc2_, eles_tight);

   edm::Handle<edm::View<pat::Electron>> eles_medium;
   iEvent.getByToken(EleSrc1_, eles_medium);
  
   edm::Handle<edm::View<pat::Electron>> eles0;
   iEvent.getByToken(EleSrc0_, eles0);
 

// ************************* Gen Level Information******************  //
	  for(size_t ik=0; ik<genParticles->size();ik++)
	{
				if(abs((*genParticles)[ik].pdgId())==24)
				{
                                   const reco::Candidate* pw = &(*genParticles)[ik];
                                     for(int i=0;pw->daughter(i)!=NULL;i++)
					{  
					const reco::Candidate* pl = pw->daughter(i);
					if( (abs(pl->pdgId())==11) || (abs(pl->pdgId())==12) || (abs(pl->pdgId())==13) || (abs(pl->pdgId())==14) ){
                                                if(abs(pl->pdgId())==11)
                                                {
                                                        gen_ele_pt=pl->pt();
                                                        gen_ele_eta=pl->eta();
                                                        gen_ele_phi=pl->phi();
                                                }
					}//end of w daugter loop
					}
				}//end of if w
         }

 	//****************  gen match, tight and medium ********************// 
        if(gen_ele_pt>0. && eles_tight->size()>0) {
            double drelematch=10000.;  size_t mk=0;
            for(size_t ik=0; ik<eles_tight->size();ik++)    {
               double drtemp=deltaR(gen_ele_eta,gen_ele_phi,(*eles_tight)[ik].eta(),(*eles_tight)[ik].phi());
               if (drtemp<drelematch) {drelematch=drtemp; mk=ik;}
            }
            gen_tight_match_ele_pt=(*eles_tight)[mk].pt();
            gen_tight_match_ele_eta=(*eles_tight)[mk].eta();
            gen_tight_match_ele_phi=(*eles_tight)[mk].phi();
            gen_tight_match_ele_dr=drelematch;
         }

        if(gen_ele_pt>0. && eles_medium->size()>0) {
            double drelematch=10000.;  size_t mk=0;
            for(size_t ik=0; ik<eles_medium->size();ik++)    {
               double drtemp=deltaR(gen_ele_eta,gen_ele_phi,(*eles_medium)[ik].eta(),(*eles_medium)[ik].phi());
               if (drtemp<drelematch) {drelematch=drtemp; mk=ik;}
            }
            gen_medium_match_ele_pt=(*eles_medium)[mk].pt();
            gen_medium_match_ele_eta=(*eles_medium)[mk].eta();
            gen_medium_match_ele_phi=(*eles_medium)[mk].phi();
            gen_medium_match_ele_dr=drelematch;
         }

	//******************** reco *********************//
       if(gen_ele_pt>0. && eles0->size()>0) {
            double drelematch=10000.;  size_t mk=0;
            for(size_t ik=0; ik<eles0->size();ik++)    {
               double drtemp=deltaR(gen_ele_eta,gen_ele_phi,(*eles0)[ik].eta(),(*eles0)[ik].phi());
               if (drtemp<drelematch) {drelematch=drtemp; mk=ik;}
            }
            reco_ele_pt=(*eles0)[mk].pt();
            reco_ele_eta=(*eles0)[mk].eta();
            reco_ele_phi=(*eles0)[mk].phi();
            reco_ele_dr=drelematch;
         }

       //********************** reco match, tight and medium****************//
       if(reco_ele_pt>0. && eles_tight->size()>0) {
            double drelematch=10000.;  size_t mk=0;
            for(size_t ik=0; ik<eles_tight->size();ik++)    {
               double drtemp=deltaR(reco_ele_eta,reco_ele_phi,(*eles_tight)[ik].eta(),(*eles_tight)[ik].phi());
               if (drtemp<drelematch) {drelematch=drtemp; mk=ik;}
            }
            reco_tight_match_ele_pt=(*eles_tight)[mk].pt();
            reco_tight_match_ele_eta=(*eles_tight)[mk].eta();
            reco_tight_match_ele_phi=(*eles_tight)[mk].phi();
            reco_tight_match_ele_dr=drelematch;
         }

       if(reco_ele_pt>0. && eles_medium->size()>0) {
            double drelematch=10000.;  size_t mk=0;
            for(size_t ik=0; ik<eles_medium->size();ik++)    {
               double drtemp=deltaR(reco_ele_eta,reco_ele_phi,(*eles_medium)[ik].eta(),(*eles_medium)[ik].phi());
               if (drtemp<drelematch) {drelematch=drtemp; mk=ik;}
            }
            reco_medium_match_ele_pt=(*eles_medium)[mk].pt();
            reco_medium_match_ele_eta=(*eles_medium)[mk].eta();
            reco_medium_match_ele_phi=(*eles_medium)[mk].phi();
            reco_medium_match_ele_dr=drelematch;
         }

       outTree_->Fill();
}
//-------------------------------------------------------------------------------------------------------------------------------------//


void EDBRTreeMaker::setDummyValues() {
     gen_ele_pt     = -99;
     gen_ele_eta    = -99;
     gen_ele_phi    = -99;
     gen_tight_match_ele_pt     = -99;
     gen_tight_match_ele_eta    = -99;
     gen_tight_match_ele_phi    = -99;
     gen_tight_match_ele_dr     =  99;
     gen_medium_match_ele_pt     = -99;
     gen_medium_match_ele_eta    = -99;
     gen_medium_match_ele_phi    = -99;
     gen_medium_match_ele_dr     =  99;
     reco_ele_pt     = -99;
     reco_ele_eta    = -99;
     reco_ele_phi    = -99;
     reco_ele_dr     =  99;
     reco_tight_match_ele_pt     = -99;
     reco_tight_match_ele_eta    = -99;
     reco_tight_match_ele_phi    = -99;
     reco_tight_match_ele_dr     =  99;
     reco_medium_match_ele_pt     = -99;
     reco_medium_match_ele_eta    = -99;
     reco_medium_match_ele_phi    = -99;
     reco_medium_match_ele_dr     =  99;

}

// ------------ method called once each job just before starting event loop  ------------
void 
EDBRTreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void EDBRTreeMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
   std::cout<<"\n*********************************************\n\n";
}

void EDBRTreeMaker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
std::cout << "EDBRTreeMaker endJob()... endRun" << std::endl;
}


void
EDBRTreeMaker::endJob() {
  std::cout << "EDBRTreeMaker endJob()..." << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(EDBRTreeMaker);
