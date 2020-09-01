/* created 2020-05-24 by Andreas Gaertner
Modified by Matthew Ens places an SDOM (modelled from an actual sDOM) 
and K40 decay electrons are simulated within a 3m radius sphere

Compile with
source geant4.sh
mkdir build
cd build
cmake ../
make

*/
// Include the various header files for the simulation
#include "FTFP_BERT.hh"
#include "G4MuonMinus.hh"
#include "G4NistManager.hh"
#include "G4OpticalPhysics.hh"
#include "G4Orb.hh"
#include "G4OpticalPhoton.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4RunManager.hh"
#include "G4SingleParticleSource.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VisExecutive.hh"

#include "G4Track.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4IonConstructor.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Transform3D.hh"
#include "G4SteppingVerbose.hh"
#define OVERLAP false
using namespace std;

/*  everything is this map can be changed using command line args
e.g. ./K40coinc WorldRadius 100 n 1
to run with GUI: ./K40coinc ui 1 */
map<string, double> var = {
	//{"DetectorRadius",   0.075},  // m no longer set since we want an accurate sDOM
	// {"DetectorDistance", 0.8},    // m no longer set since we want an accurate sDOM
	{"WorldRadius",      20.},      // m
	{"RefractiveIndex",  1.42},
	{"AbsorptionLength", 30},     // m
	{"ScatteringLength", 60},     // m
	{"n",                0},      // number of events
	{"ui",               0},
	{"Seed",             0},
	{"Verbose",          0},
	{"angle",            0}
};
map<string, G4ThreeVector> volumePositions = {};

//This is the detector definition class, it defines the world full of water along with the simulated sDOM in the centre of the world

class K40coincDetectorConstruction : public G4VUserDetectorConstruction{
public:
	K40coincDetectorConstruction(): G4VUserDetectorConstruction(){}
    virtual G4VPhysicalVolume* Construct(){

		// water properties
		G4NistManager* nist = G4NistManager::Instance();
		G4Material* water = nist->FindOrBuildMaterial("G4_WATER");
		const size_t entries = 2;
		G4double PhotonEnergy     [entries] = {1*eV, 4.5*eV};
		G4double RefractiveIndex  [entries] = {var["RefractiveIndex"], var["RefractiveIndex"]};
		//G4double AbsorptionLength [entries] = {var["AbsorptionLength"]*m, var["AbsorptionLength"]*m};
		G4double ScatteringLength [entries] = {var["ScatteringLength"]*m, var["ScatteringLength"]*m};

		const size_t absEntries = 8;
		G4double absEnergies      [absEntries] = {1*eV, 2.32834*eV, 2.32835*eV, 2.85021*eV, 2.85022*eV, 3.22037*eV, 3.22038*eV,4.5*eV};
		G4double AbsorptionLength [absEntries] = {5.*m, 5.*m,31.87*m,31.87*m,17.56*m,17.56*m,9.21*m,9.21*m};

		G4MaterialPropertiesTable* wp = new G4MaterialPropertiesTable();
		wp->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex, entries);
		wp->AddProperty("ABSLENGTH", absEnergies, AbsorptionLength, absEntries);
		wp->AddProperty("RAYLEIGH", PhotonEnergy, ScatteringLength, entries);
		water->SetMaterialPropertiesTable(wp);

		// geometry
		G4double worldRadius = var["WorldRadius"]*m;
		G4Orb* world = new G4Orb("World", worldRadius);
		G4LogicalVolume* world_logic = new G4LogicalVolume(world, water, "World");
		G4VPhysicalVolume* world_phys = new G4PVPlacement(0, G4ThreeVector(), world_logic, 
		                                                  "World", 0, false, 0, true);

		//make the titanium housing
		G4Material* titanium = nist->FindOrBuildMaterial("G4_Ti");
		G4ThreeVector centre = G4ThreeVector(0,0,0);
		G4String tubename = "housing"; //below: (name,r_in,r_out,0.5*length,start_ang,end_ang)
		G4Tubs* housing = new G4Tubs(tubename, 0.,0.065*m,0.2441*m,0.*degree,360*degree);
		G4LogicalVolume* housing_logic = new G4LogicalVolume(housing,titanium,tubename);
		volumePositions[tubename] = centre;
		new G4PVPlacement(0,centre,housing_logic,tubename,world_logic, false, 0, OVERLAP);

		G4ThreeVector x = G4ThreeVector(0, 0, 0.2441*m);
		G4RotationMatrix rot = G4RotationMatrix();
		rot.rotateX(180*degree);
		G4Transform3D transform = G4Transform3D(rot,-x);

		
		//make outer glass spheres
		G4Material* glass = nist->FindOrBuildMaterial("G4_GLASS_PLATE");
		G4double glassRefractiveIndex [entries] = {1.52,1.52};
		G4MaterialPropertiesTable* glassp = new G4MaterialPropertiesTable();
		glassp->AddProperty("RINDEX", PhotonEnergy, glassRefractiveIndex, entries);
		glass->SetMaterialPropertiesTable(glassp);

		G4String lGlassSphere_name = "lSphere";
		G4Sphere* lSphere = new G4Sphere(lGlassSphere_name,0.0381*m,0.05715*m,0*degree,360*degree,0*degree,90*degree);
		G4LogicalVolume* lSphere_logic = new G4LogicalVolume(lSphere,glass,lGlassSphere_name);
		volumePositions[lGlassSphere_name] = x;
		new G4PVPlacement(0,x,lSphere_logic, lGlassSphere_name, world_logic, false, 0, OVERLAP);

		G4String rGlassSphere_name = "rSphere";
		G4Sphere* rSphere = new G4Sphere(rGlassSphere_name,0.0381*m,0.05715*m,0*degree,360*degree,0*degree,90*degree);
		G4LogicalVolume* rSphere_logic = new G4LogicalVolume(rSphere,glass,rGlassSphere_name);
		volumePositions[rGlassSphere_name] = -x;
		new G4PVPlacement(transform, rSphere_logic, rGlassSphere_name, world_logic, false, 0, OVERLAP);

		G4double detPhi = 53.13; //actual is 53.13

		G4String eHousel_name = "housing";
		G4Sphere* eHousel = new G4Sphere(eHousel_name,0.,0.0381*m,0.,360*degree,detPhi*deg,(90-detPhi)*deg);
		G4LogicalVolume* eHousel_logic = new G4LogicalVolume(eHousel, titanium, eHousel_name);
		volumePositions[eHousel_name] = x;
		new G4PVPlacement(0, x, eHousel_logic, eHousel_name, world_logic, false, 0, OVERLAP);

		G4String eHouser_name = "housing";
		G4Sphere* eHouser = new G4Sphere(eHouser_name,0*m,0.0381*m,0.,360*degree,detPhi*deg,(90-detPhi)*deg);
		G4LogicalVolume* eHouser_logic = new G4LogicalVolume(eHouser, titanium, eHouser_name);
		volumePositions[eHouser_name] = -x;
		new G4PVPlacement(transform, eHouser_logic, eHouser_name, world_logic, false, 0, OVERLAP);

		//make pmts
		G4String lPMT_name = "up"; //(name,r_in,r_out,phi_min,phi_max,th_min,th_max)
		G4Sphere* lPMT = new G4Sphere(lPMT_name, 0.*m,0.0381*m,0*degree,360*degree,0,detPhi*deg);
		G4LogicalVolume* lPMT_logic = new G4LogicalVolume(lPMT, glass, lPMT_name);
		volumePositions[lPMT_name] = x;
		new G4PVPlacement(0, x, lPMT_logic, lPMT_name, world_logic, false, 0, OVERLAP);
		
		G4String rPMT_name = "down";
		G4Sphere* rPMT = new G4Sphere(rPMT_name, 0.*m,0.0381*m,0*degree,360*degree,0,detPhi*deg);
		G4LogicalVolume* rPMT_logic = new G4LogicalVolume(rPMT, glass, rPMT_name);
		volumePositions[rPMT_name] = -x;
		new G4PVPlacement(transform, rPMT_logic, rPMT_name, world_logic, false, 0, OVERLAP);

		return world_phys;
	}
};

//This is the class for the K40 beta decay electrons it places them randomly in the world with the proper energy
class K40coincPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{
public:
	K40coincPrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(){
		// electron source from 40K decay
		//THIS WAS COMMENTED OUT!
		G4ParticleDefinition* part = G4ParticleTable::GetParticleTable()->FindParticle( "e-" );

		genSource = new G4GeneralParticleSource();
		G4SingleParticleSource* source = genSource->GetCurrentSource();
		source->SetParticleDefinition(part);
		G4SPSPosDistribution* position = source->GetPosDist();
		position->SetPosDisType("Point");

		G4cout << "HIHO" << G4RandFlat::shoot() << G4endl;
		G4SPSAngDistribution* angle = source->GetAngDist();
		angle->SetAngDistType("iso");
		angle->DefineAngRefAxes("angref1", G4ThreeVector(1., 0., 0.));
		angle->DefineAngRefAxes("angref2", G4ThreeVector(0., 1., 0.));
		angle->SetMinTheta(0.*deg);
		angle->SetMaxTheta(180.*deg);
		angle->SetMinPhi(0.*degree);
		angle->SetMaxPhi(360.*degree);
		G4SPSEneDistribution* energy = source->GetEneDist();
		energy->SetEnergyDisType("Gauss");
		energy->SetMonoEnergy(0.7*MeV); 
		energy->SetEmin(0.05*MeV);
		energy->SetEmax(1.31*MeV);
		energy->SetBeamSigmaInE(0.4*MeV);

	}
	virtual ~K40coincPrimaryGeneratorAction(){
		delete genSource;
	}
	virtual void GeneratePrimaries(G4Event* event){
		
		//If want entire world
		G4double u = G4RandGauss::shoot(0,1);
		G4double v = G4RandGauss::shoot(0,1);
		G4double w = G4RandGauss::shoot(0,1);
		G4double r = var["WorldRadius"]*std::pow(G4UniformRand(),1./3);

		G4double normalization = std::pow((u*u + v*v + w*w),0.5);

		G4ThreeVector randpos = r*(G4ThreeVector(u,v,w))/normalization*m;
		G4SingleParticleSource* source = genSource->GetCurrentSource();
    	source->SetParticlePolarization(G4RandomDirection());
    	source->SetParticleTime(G4UniformRand()*ms);
    	source->GetPosDist()->SetCentreCoords(randpos);
    	source->GeneratePrimaryVertex(event);

	}
private:
	G4GeneralParticleSource* genSource;
};

//This class adds more verbose output if the verbose option is chosen
class K40coincSteppingVerbose : public G4SteppingVerbose {
public:
	K40coincSteppingVerbose() : G4SteppingVerbose(){}
	virtual void TrackingStarted(){
		CopyState();
  
	  G4int prec = G4cout.precision(3);
	  
	  //Step zero
	  //  
	  if( verboseLevel > 0 ){
		   	G4cout << std::setw( 5) << "Step#"      << " "
	           << std::setw( 6) << "X"          << "    "
	           << std::setw( 6) << "Y"          << "    "  
	           << std::setw( 6) << "Z"          << "    "
	           << std::setw( 9) << "KineE"      << " "
	           << std::setw( 9) << "dEStep"     << " "  
	           << std::setw( 9) << "Time"		<< " "
	           << std::setw(10) << "StepLeng"  
	           << std::setw(10) << "TrakLeng"
	           << std::setw(10) << "Volume"     << "  "
	           << std::setw(10) << "Process"    << G4endl;             

		    G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
		        << std::setw(6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
		        << std::setw(6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
		        << std::setw(6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
		        << std::setw(6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
		        << std::setw(6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
		        << std::setw(6) << G4BestUnit(fTrack->GetGlobalTime(), "Time")
		        << std::setw(6) << G4BestUnit(fStep->GetStepLength(),"Length")
		        << std::setw(6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
		        << std::setw(10) << fTrack->GetVolume()->GetName()
		        << "   initStep" << G4endl;        
	  	}
	  G4cout.precision(prec);
	}	

	virtual void StepInfo(){

		CopyState();
    
	  G4int prec = G4cout.precision(3);

	  if( verboseLevel >= 1 ){
	    	if( verboseLevel >= 4 ) VerboseTrack();
	    	if( verboseLevel >= 3 ){
		      G4cout << G4endl;    
		      G4cout << std::setw( 5) << "#Step#"     << " "
	             << std::setw( 6) << "X"          << "    "
	             << std::setw( 6) << "Y"          << "    "  
	             << std::setw( 6) << "Z"          << "    "
	             << std::setw( 9) << "KineE"      << " "
	             << std::setw( 9) << "dEStep"     << " "  
	             << std::setw( 9) << "Time"       << " "
	             << std::setw(10) << "StepLeng"     
	             << std::setw(10) << "TrakLeng" 
	             << std::setw(10) << "Volume"    << "  "
	             << std::setw(10) << "Process"   << G4endl;                  
	    	}

		    G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
		        << std::setw(6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
		        << std::setw(6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
		        << std::setw(6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
		        << std::setw(6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
		        << std::setw(6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
		        << std::setw(6) << G4BestUnit(fTrack->GetGlobalTime(), "Time")
		        << std::setw(6) << G4BestUnit(fStep->GetStepLength(),"Length")
		        << std::setw(6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
		        << std::setw(10) << fTrack->GetVolume()->GetName();


		    const G4VProcess* process 
		                      = fStep->GetPostStepPoint()->GetProcessDefinedStep();
		    G4String procName = " UserLimit";
		    if (process) procName = process->GetProcessName();
		    if (fStepStatus == fWorldBoundary) procName = "OutOfWorld";
		    G4cout << "   " << std::setw(10) << procName;
		    G4cout << G4endl;
		}
	    G4cout.precision(prec);
	}
};


//This class is the detector response, if you want extra data on a hit to be output you would add it here
class K40coincSteppingAction : public G4UserSteppingAction{
public:
	K40coincSteppingAction() : G4UserSteppingAction(){}
	virtual void UserSteppingAction(const G4Step* step){
		// detect photons, define variables
		G4StepPoint*      startpoint = step->GetPreStepPoint();
		G4TouchableHandle geometry   = startpoint->GetTouchableHandle();
		G4String          volumeName = geometry->GetVolume()->GetName();
		static bool       dataStart  = true;
		static bool       first      = true;
		static int        up         = 0, down = 0; 
		static G4double   quantum_eff = 1., uptime, downtime, uplen, downlen,dene1,uene1;

		bool record_hit = (G4UniformRand()<quantum_eff);

		/*if (volumeName=="World" && step->GetTrack()->GetParticleDefinition()->GetParticleName()=="opticalphoton"){
			
			G4double wene = step->GetTrack()->GetKineticEnergy()/eV;
			
			if (wene > 5. || wene < 1.){
				step->GetTrack()->SetTrackStatus(fStopAndKill);
			}
		}*/

		//record a hit if we a photon is found in the 'up' or 'down' volume
		if (volumePositions.count(volumeName) > 0){ 
			if (record_hit && (volumeName=="up" || volumeName =="down")){
				if (dataStart){
					G4cout << "START" << G4endl; 
					dataStart=false; 
					G4cout << std::fixed;
					G4cout << std::setprecision(1);
				}
				//if (first){ G4cout << "#DET\tdetector\ttime[ns]" << G4endl; first = false;}
				//G4double time = startpoint->GetGlobalTime()/ns;
				if (volumeName == "up" && up++ == 0){
					uptime = startpoint->GetGlobalTime()/ns;
					G4ThreeVector ucreatedPos = step->GetTrack()->GetVertexPosition();
					uplen = ucreatedPos.mag()/cm;
					G4double uene = step->GetTrack()->GetKineticEnergy()/eV;
					if (up>1){G4cout << "#\t1\t0\t" << uptime << "\t" << uplen<< "\t" << uene<<G4endl;}
					else {uene1 = uene;}				
				}
				if (volumeName == "down" && down++ == 0){
					downtime = startpoint->GetGlobalTime()/ns;
					G4ThreeVector dcreatedPos = step->GetTrack()->GetVertexPosition();
					downlen = dcreatedPos.mag()/cm;
					G4double dene = step->GetTrack()->GetKineticEnergy()/eV;
					if(down>1){G4cout << "#\t5\t0\t" << downtime << "\t" << downlen<< "\t" << dene<<G4endl;}
					else{dene1=dene;}
				}
				//G4ThreeVector createdPos = step->GetTrack()->GetVertexPosition();
				//G4double tracklen = createdPos.mag();
				//G4cout << volumeName << "\t" << time   << G4endl;
				first=false;
			}
			if (volumeName=="up" || volumeName == "down" || volumeName == "housing"){
				step->GetTrack()->SetTrackStatus(fStopAndKill);
			}
		}
//output number of photons detected once we start the next electron simulation, reset all counts
		if (step->GetTrack()->GetParentID()==0 && step->GetTrack()->GetCurrentStepNumber()==1 && !(first)){
			if (up>0){G4cout << "1" << "\t" << "0" << "\t" << uptime << "\t" << up << "\t" << uplen << "\t" << uene1<< G4endl;}
			if (down>0){G4cout << "5" << "\t" << "0" << "\t" << downtime << "\t" << down << "\t" << downlen <<"\t"<<dene1<<G4endl;}
			//if (up>0 && down>0) {G4cout << "coincidence!" << G4endl;}
			//G4cout << "\nNew Track!" << G4endl;
			first = true;
			up = 0; down = 0;

		}
	}

};
//main function which gets the physics together and runs everything
int main(int argc,char** argv){
	for (int i = 1; i < argc; i += 2){
		if (var.find(string(argv[i])) == var.end()){
			throw invalid_argument("Error parsing args. Unknown parameter.");}
		try {
			var[string(argv[i])] = stod(argv[i+1]);
		}catch (const std::exception& e){
			throw invalid_argument("Error parsing args. Value missing.");}}
	for (map<string, double>::iterator it = var.begin(); it != var.end(); it++)
		G4cout << "VAR " << it->first << " " << it->second << G4endl;

	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(var["Seed"]);
	G4RandFlat::setTheEngine(new CLHEP::RanecuEngine);
	G4RandFlat::setTheSeed(var["Seed"]);
	G4RandGauss::setTheEngine(new CLHEP::RanecuEngine);
	G4RandGauss::setTheSeed(var["Seed"]);
	G4RandExponential::setTheEngine(new CLHEP::RanecuEngine);
	G4RandExponential::setTheSeed(var["Seed"]);
	CLHEP::RandPoisson::setTheEngine(new CLHEP::RanecuEngine);
	CLHEP::RandPoisson::setTheSeed(var["Seed"]);

	G4VSteppingVerbose::SetInstance(new K40coincSteppingVerbose);


	G4RunManager* runManager = new G4RunManager;
	runManager->SetUserInitialization(new K40coincDetectorConstruction());
	G4VModularPhysicsList* p = new FTFP_BERT;
	G4OpticalPhysics* o = new G4OpticalPhysics();
	o->SetMaxNumPhotonsPerStep(100);
	o->SetMaxBetaChangePerStep(0.1);
	o->SetTrackSecondariesFirst(kCerenkov, true);
	p->RegisterPhysics(o);

	runManager->SetUserInitialization(p);
	runManager->SetUserAction(new K40coincPrimaryGeneratorAction());
	runManager->SetUserAction(new K40coincSteppingAction());

	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();
	G4UImanager* u = G4UImanager::GetUIpointer();
	u->ApplyCommand("/run/initialize");

	if (var["ui"] != 0.0){
		G4UIExecutive* ui = new G4UIExecutive(argc, argv);
		u->ApplyCommand("/vis/open OGL 600x600-0+0");
		u->ApplyCommand("/vis/drawVolume");
		u->ApplyCommand("/vis/viewer/set/lineSegmentsPerCircle 15");
		u->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 90 0");
		u->ApplyCommand("/vis/scene/add/trajectories smooth");
		u->ApplyCommand("/vis/scene/endOfEventAction accumulate");
		u->ApplyCommand("/vis/viewer/set/style wireframe");
		//u->ApplyCommand("/tracking/verbose 1");
		u->ApplyCommand(G4String("/run/beamOn ") +G4String(to_string((int)var["n"])));
		ui->SessionStart();
		delete ui;
	}
	else if (var["Verbose"] != 0.0){
		u->ApplyCommand("/tracking/verbose 1");
		u->ApplyCommand(G4String("/run/beamOn ") +G4String(to_string((int)var["n"])));
	}
	else{
		//min=3867 s-1m-3; max=8801 s-1m-3 -> 4.93414*RAND+3.867465 ms-1 m-3
		//actual: 12 133 s-1m-3 -> 12.00*RAND+12.133 ms-1 m-3
		G4double worldVol = 4/3*pi*std::pow(var["WorldRadius"],3);
		G4cout << "World Volume: " << worldVol << G4endl;
		//G4double randnum = 12.00*worldVol*G4UniformRand()+12.133*worldVol; /*500.*G4UniformRand()+144.;*/ 
		G4double pInVol = 12.133*worldVol; //ms-1
		G4long numparticles = CLHEP::RandPoisson::shoot(pInVol);
		//G4double randnum = randRate*worldVol;
		//G4int numparticles = randnum;
		G4cout << "Particles generated: " << numparticles << G4endl;
		u->ApplyCommand(G4String("/run/beamOn ") +G4String(to_string((int)numparticles)));
	}
	delete visManager;
	delete runManager;
}
