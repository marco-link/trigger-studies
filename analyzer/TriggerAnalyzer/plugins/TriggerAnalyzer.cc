// -*- C++ -*-
//
// Package:    analyzer/TriggerAnalyzer
// Class:      TriggerAnalyzer
// 
/**\class TriggerAnalyzer TriggerAnalyzer.cc analyzer/TriggerAnalyzer/plugins/TriggerAnalyzer.cc

 Description: analyzer to extract data for trigger efficiency plots for the B2G dijet full hadronic analysis

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marco Link
//         Created:  Wed, 20 Sep 2017 07:37:29 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/Jet.h"



class TriggerAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TriggerAnalyzer(const edm::ParameterSet&);
      ~TriggerAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        bool passedPreselections(const edm::Event& iEvent);
        bool tightJetID(const pat::Jet& j);


        // parameter
        bool DoPreselection;
        bool DoFilter;
        bool DoTightJetID;

        double EtaCut;
        double PtCut;
        double dEtaCut;
        double HTEtaCut;
        double HTPtCut;


        // input
        edm::EDGetTokenT<reco::VertexCollection>    vtxToken_       ;
        edm::Handle<reco::VertexCollection>         vertices_       ;

        edm::EDGetTokenT<pat::JetCollection>        jetInputToken_  ;
        edm::Handle<pat::JetCollection>             jets_           ;

        edm::EDGetTokenT<pat::JetCollection>        htjetInputToken_;
        edm::Handle<pat::JetCollection>             htjets_         ;


        // trigger stuff
        edm::EDGetTokenT<edm::TriggerResults>                       triggerToken_    ;
        edm::Handle< edm::TriggerResults>                           HLTtriggers_     ;

        edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>    triggerObjects_  ;
        edm::Handle<pat::TriggerObjectStandAloneCollection>         triggerObjects   ;

        edm::EDGetTokenT<edm::TriggerResults>                       noiseFilterToken_;
        edm::Handle< edm::TriggerResults>                           noiseFilterBits_ ;



        edm::EDGetTokenT<bool> HBHENoiseFilterLoose_Selector_   ;
        edm::EDGetTokenT<bool> HBHENoiseFilterTight_Selector_   ;
        edm::EDGetTokenT<bool> HBHENoiseIsoFilter_Selector_     ;

        std::string ECALDeadCellNoiseFilter_Selector_               ;
        std::string EEBadScNoiseFilter_Selector_                    ;
        std::string globalTightHalo2016Filter_Selector_             ;
        std::string muonBadTrackFilter_Selector_                    ;
        std::string chargedHadronTrackResolutionFilter_Selector_    ;

        unsigned int passed_filter  ;
        unsigned int found_jets     ;
        unsigned int passed_etacut  ;

        pat::Jet leadingJet     ;
        pat::Jet subleadingJet  ;

        std::vector<std::string> triggerlist;
        std::ofstream out;
};


TriggerAnalyzer::TriggerAnalyzer(const edm::ParameterSet& ps)
{
    //now do what ever initialization is needed
    usesResource("TFileService");
    // Get parameters from configuration file
    vtxToken_ = consumes<reco::VertexCollection>(ps.getParameter<edm::InputTag>("vertices"));


    jetInputToken_      = consumes<pat::JetCollection>(ps.getParameter<edm::InputTag>("jets"));
    htjetInputToken_      = consumes<pat::JetCollection>(ps.getParameter<edm::InputTag>("htjets"));
    noiseFilterToken_   = consumes<edm::TriggerResults>(ps.getParameter<edm::InputTag>("NoiseFilter"));;
    triggerToken_       = consumes<edm::TriggerResults>(ps.getParameter<edm::InputTag>("TriggerResults"));
    triggerObjects_     = consumes<pat::TriggerObjectStandAloneCollection>(ps.getParameter<edm::InputTag>("TriggerObjects"));

    HBHENoiseFilterLoose_Selector_   = consumes<bool>(ps.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseFilterLoose"))  ;
    HBHENoiseFilterTight_Selector_   = consumes<bool>(ps.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseFilterTight"))  ;
    HBHENoiseIsoFilter_Selector_     = consumes<bool>(ps.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseIsoFilter"))    ;

    ECALDeadCellNoiseFilter_Selector_           = ps.getParameter<std::string> ("noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter") ;
    EEBadScNoiseFilter_Selector_                = ps.getParameter<std::string> ("noiseFilterSelection_eeBadScFilter")                      ;
    globalTightHalo2016Filter_Selector_         = ps.getParameter<std::string> ("noiseFilterSelection_globalTightHalo2016Filter")          ;
    muonBadTrackFilter_Selector_                = ps.getParameter<std::string> ("noiseFilterSelection_muonBadTrackFilter")                 ;
    chargedHadronTrackResolutionFilter_Selector_= ps.getParameter<std::string> ("noiseFilterSelection_chargedHadronTrackResolutionFilter") ;

    // parameter
    DoPreselection  = ps.getUntrackedParameter<bool> ("DoPreselection");
    DoFilter        = ps.getUntrackedParameter<bool> ("DoFilter");
    DoTightJetID    = ps.getUntrackedParameter<bool> ("DoTightJetID");

    EtaCut          = ps.getUntrackedParameter<double>("EtaCut");
    PtCut           = ps.getUntrackedParameter<double>("PtCut");
    dEtaCut         = ps.getUntrackedParameter<double>("dEtaCut");
    HTEtaCut        = ps.getUntrackedParameter<double>("HTEtaCut");
    HTPtCut         = ps.getUntrackedParameter<double>("HTPtCut");

    passed_filter   = 0;
    found_jets      = 0;
    passed_etacut   = 0;


    std::istringstream names(ps.getUntrackedParameter<std::string>("triggernames"));
    std::string token;

    while (std::getline(names, token, ','))
    {
        triggerlist.push_back(token);
    }

    if(triggerlist.size()<2)
    {
        throw cms::Exception("triggerlist not found") << "need at least 2 trigger for trigger studies (nominator/denominator), found only " << triggerlist.size() << std::endl;
    }

    out.open(ps.getUntrackedParameter<std::string>("target"), std::ios::out);

    if(!out.is_open())
    {
        throw cms::Exception("outfile not generated") << "cant open outfile: " << ps.getUntrackedParameter<std::string>("target") << std::endl;
    }

    out << "run,ID,lumi,Mjj,HT,pt1,pt2,phi1,phi2,eta1,eta2,softdrop1,softdrop2";

    for(unsigned int i = 0; i<triggerlist.size(); i++)
    {
        out << "," << triggerlist[i];
    }
    out << std::endl;
}


TriggerAnalyzer::~TriggerAnalyzer()
{
   out.close();
   std::cout << std::endl << "- preselections -" << std::endl << "passed filter: " << passed_filter << "\nfound good jets: " << found_jets << "\npassed deltaEta-cut: " << passed_etacut << std::endl;
}


void TriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    iEvent.getByToken(vtxToken_, vertices_);
    iEvent.getByToken(jetInputToken_, jets_);


    if(!DoPreselection || passedPreselections(iEvent))
    {
        iEvent.getByToken(triggerToken_, HLTtriggers_);
        iEvent.getByToken(triggerObjects_  , triggerObjects);

        const edm::TriggerNames& trigNames = iEvent.triggerNames(*HLTtriggers_);
        std::vector<unsigned int> fired;
        for(unsigned int i = 0; i<triggerlist.size(); i++)
        {
            fired.push_back(0);
        }


        bool Any = false;

        // loop over trigger
        for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {
            if(HLTtriggers_->accept(i))
            {
                for(unsigned int j = 0; j<triggerlist.size(); j++)
                {
                    if(trigNames.triggerName(i).find(triggerlist[j]) != std::string::npos)
                    {
                        fired[j] = 1;
                        Any = true;
                    }
                }
            }
        }

        if(Any)
        {
            double Mjj = (leadingJet.p4() + subleadingJet.p4()).M();
            double HT = 0;
            iEvent.getByToken(htjetInputToken_, htjets_);
            for (auto htjet = htjets_->begin(); htjet != htjets_->end(); htjet++)
            {
                if(htjet->p4().pt() > HTPtCut && fabs(htjet->p4().Eta()) < HTEtaCut)
                {
                    HT = HT + htjet->p4().pt();
                }
            }

            double pt1 = leadingJet.p4().pt();
            double pt2 = subleadingJet.p4().pt();

            double phi1 = leadingJet.p4().Phi();
            double phi2 = subleadingJet.p4().Phi();

            double eta1 = leadingJet.p4().Eta();
            double eta2 = subleadingJet.p4().Eta();

            double softdrop1 = leadingJet.userFloat("ak8PFJetsPuppiSoftDropMass");
            double softdrop2 = subleadingJet.userFloat("ak8PFJetsPuppiSoftDropMass");

            out << iEvent.id().run()
                << ", " << iEvent.id().event()
                << ", " << iEvent.id().luminosityBlock()
                << ", " << Mjj
                << ", " << HT
                << ", " << pt1
                << ", " << pt2
                << ", " << phi1
                << ", " << phi2
                << ", " << eta1
                << ", " << eta2
                << ", " << softdrop1
                << ", " << softdrop2;
            for(unsigned int i = 0; i < fired.size(); i++)
            {
                out << ", " << fired.at(i);
            }
            out << std::endl;
        }
    }
}


void TriggerAnalyzer::beginJob()
{
}


void TriggerAnalyzer::endJob() 
{
}


void
TriggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


bool sortByPt(const pat::Jet &lhs, const pat::Jet &rhs){return lhs.p4().pt() > rhs.p4().pt();};
bool TriggerAnalyzer::passedPreselections(const edm::Event& iEvent)
{
    // primary vertex
    reco::VertexCollection::const_iterator firstGoodVertex = vertices_->end();
    int firstGoodVertexIdx = 0;
    for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx, ++firstGoodVertexIdx)
    {
        bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
        if( !isFake && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0)
        {
            firstGoodVertex = vtx;
            break;
        }
    }

    if (firstGoodVertex==vertices_->end())
    {
        return false;
    }


    bool passFilter = true;

    edm::Handle<bool> HBHENoiseIsoFilterResultHandle;
    iEvent.getByToken(HBHENoiseIsoFilter_Selector_, HBHENoiseIsoFilterResultHandle);
    passFilter = *HBHENoiseIsoFilterResultHandle;
    if (!HBHENoiseIsoFilterResultHandle.isValid())
    {
        LogDebug("") << "CaloTowerAnalyzer: Could not find HBHENoiseFilterResult" << std::endl;
    }

    edm::Handle<bool> HBHENoiseFilterLooseResultHandle;
    iEvent.getByToken(HBHENoiseFilterLoose_Selector_, HBHENoiseFilterLooseResultHandle);
    passFilter = passFilter && *HBHENoiseFilterLooseResultHandle;
    if (!HBHENoiseFilterLooseResultHandle.isValid())
    {
        LogDebug("") << "CaloTowerAnalyzer: Could not find HBHENoiseFilterResult" << std::endl;
    }


    iEvent.getByToken(noiseFilterToken_, noiseFilterBits_);
    const edm::TriggerNames &names = iEvent.triggerNames(*noiseFilterBits_);


    for (unsigned int i = 0, n = noiseFilterBits_->size(); i < n; ++i)
    {
        if (names.triggerName(i) == ECALDeadCellNoiseFilter_Selector_)
        {
            passFilter = passFilter && noiseFilterBits_->accept(i); // under scrutiny
        }

        if (names.triggerName(i) == EEBadScNoiseFilter_Selector_)
        {
            passFilter = passFilter && noiseFilterBits_->accept(i); // under scrutiny
        }

        if (names.triggerName(i) == globalTightHalo2016Filter_Selector_           )
        {
            passFilter = passFilter && noiseFilterBits_->accept(i); // TO BE USED FOR ICHEP 2016
        }

        if (names.triggerName(i) == muonBadTrackFilter_Selector_                  )
        {
            passFilter = passFilter && noiseFilterBits_->accept(i); // TO BE USED FOR ICHEP 2016
        }

        if (names.triggerName(i) == chargedHadronTrackResolutionFilter_Selector_  )
        {
            passFilter = passFilter && noiseFilterBits_->accept(i); // TO BE USED FOR ICHEP 2016
        }
    }

    if(!(passFilter) && DoFilter)
    {
        return false;
    }

    passed_filter++;

    pat::JetCollection Vcand;
    // filter jets
    for (auto jet = jets_->begin(); jet != jets_->end(); jet++)
    {
        if(jet->p4().pt() > PtCut && fabs(jet->p4().Eta()) < EtaCut )
        {
            Vcand.push_back(*jet);
        }
    }

    // select leading jets
    if(Vcand.size() < 2)
    {
        return false;
    }

    // sort jets (should already be sorted, just to be sure)
    sort(Vcand.begin(), Vcand.end(), sortByPt);

    leadingJet = Vcand.at(0);
    subleadingJet = Vcand.at(1);

    // tight jetID
    if(DoTightJetID && !(tightJetID(leadingJet) && tightJetID(subleadingJet)))
    {
        return false;
    }


    found_jets++;

    // deltaEta cut
    if(fabs(leadingJet.p4().Eta() - subleadingJet.p4().Eta()) > dEtaCut)
    {
        return false;
    }

    passed_etacut++;

    return true;
}


bool TriggerAnalyzer::tightJetID(const pat::Jet& j)
{
    //In sync with: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    double eta = j.eta();		
    double chf = j.chargedHadronEnergyFraction();
    double nhf = j.neutralHadronEnergyFraction(); // + j.HFHadronEnergyFraction();
    double muf = j.muonEnergy()/(j.jecFactor(0) * j.energy());  
    double nemf = j.neutralEmEnergyFraction();
    double cemf = j.chargedEmEnergyFraction();
    int chMult = j.chargedMultiplicity();
    int neMult = j.neutralMultiplicity();
    int npr    = chMult + neMult;
    int NumConst = npr;

    if(fabs(eta) <= 2.7)
    {
        return (nhf<0.90 && nemf<0.90 && NumConst>1 && muf<0.8) && ((fabs(eta)<=2.4 && chf>0 && chMult>0 && cemf<0.90) || fabs(eta)>2.4);  		
    }
    else if(fabs(eta) <= 3.0)
    {
        return (nhf<0.98 && nemf>0.01 && neMult>2);
    }
    else
    {
        return (nemf<0.90 && neMult>10);
    }
}



//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzer);