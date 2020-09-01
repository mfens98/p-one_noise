/* created 2020-05-24 by Andreas Gaertner
places an idealised SDOM (two spheres (PMTs) with isotropic
angular response) and idealised light source (isotropic)

Compile with
source geant4.sh
mkdir build
cd build
cmake ../
make

*/

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

#include "G4ParticleTable.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#define OVERLAP false
using namespace std;

/*  everything is this map can be changed using command line args
e.g. ./K40coinc WorldRadius 100 n 1
to run with GUI: ./K40coinc ui 1 */
map<string, double> var = {
	//{"DetectorRadius",   0.075},  // m no longer set since we want an accurate sDOM
	// {"DetectorDistance", 0.8},    // m no longer set since we want an accurate sDOM
	{"WorldRadius",      3},      // m
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

class K40coincDetectorConstruction : public G4VUserDetectorConstruction{
public:
	K40coincDetectorConstruction(): G4VUserDetectorConstruction(){}
    virtual G4VPhysicalVolume* Construct(){

		// water properties
		G4NistManager* nist = G4NistManager::Instance();
		G4Material* water = nist->FindOrBuildMaterial("G4_WATER");
		const size_t entries = 2;
		G4double PhotonEnergy     [entries] = {1*eV, 10*eV};
		G4double RefractiveIndex  [entries] = {var["RefractiveIndex"], var["RefractiveIndex"]};
		G4double AbsorptionLength [entries] = {var["AbsorptionLength"]*m, var["AbsorptionLength"]*m};
		G4double ScatteringLength [entries] = {var["ScatteringLength"]*m, var["ScatteringLength"]*m};
		G4MaterialPropertiesTable* wp = new G4MaterialPropertiesTable();
		wp->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex, entries);
		wp->AddProperty("ABSLENGTH", PhotonEnergy, AbsorptionLength, entries);
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

		G4double detPhi = 53.13; //actual is 53.13deg

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

class K40coincPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{
public:
	K40coincPrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(){
		// electron source from 40K decay
		
		//G4ParticleDefinition* part = G4ParticleTable::GetParticleTable()->FindParticle( "G4OpticalPhoton" );

		source = new G4SingleParticleSource();
		source->SetParticleDefinition(G4OpticalPhoton::Definition()); //source->SetParticleDefinition(part);
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
		energy->SetMonoEnergy(2.66*eV);
		energy->SetBeamSigmaInE(7.56e-2*eV);
		//energy->SetMonoEnergy(0.7*MeV); 
		//energy->SetEmin(0.05*MeV); // 465 nm
		//energy->SetEmax(1.31*MeV);
		//energy->SetBeamSigmaInE(0.4*MeV);
	}
	virtual ~K40coincPrimaryGeneratorAction(){
		delete source;
	}
	virtual void GeneratePrimaries(G4Event* event){
		// pulse shape
		G4double FWHM = 0.01*ns;
		//if want fixed postition
		G4ThreeVector detpos = G4ThreeVector(0,0,0.25*m);
		G4double fixedtheta = var["angle"]*deg;
		G4ThreeVector randpos = G4ThreeVector(0,0.25*std::sin(fixedtheta)*m,0.25*std::cos(fixedtheta)*m);
		/*//If want entire world
		G4double randR = 3*G4UniformRand();
		G4double randphi = (360*G4UniformRand())*degree;
		G4double randtheta = 180*G4UniformRand()*degree; 

		G4double randx = randR*std::sin(randtheta)*std::cos(randphi);
		G4double randy = randx*std::tan(randphi);
		G4double randz = randR*std::cos(randtheta);

		G4ThreeVector randpos = G4ThreeVector(randx*m,randy*m,randz*m);*/

		/* //If want middle section
		G4double randR = 2*std::sqrt(2)*G4UniformRand();
		G4double randtheta = 180*G4UniformRand()*degree;
		G4double randz = 2*G4UniformRand()-1;

		G4double randx = randR*std::cos(randtheta);
		G4double randy = randR*std::sin(randtheta);

		G4ThreeVector randpos = G4ThreeVector(randx*m,randy*m,randz*m);*/

    	source->SetParticlePolarization(G4RandomDirection());
    	source->SetParticleTime(G4RandGauss::shoot(0, FWHM/2.355));
    	source->GetPosDist()->SetCentreCoords(randpos+detpos);
    	source->GeneratePrimaryVertex(event);
	}
private:
	G4SingleParticleSource* source;
};

class K40coincSteppingAction : public G4UserSteppingAction{
public:
	K40coincSteppingAction() : G4UserSteppingAction(){}
	virtual void UserSteppingAction(const G4Step* step){
		// detect photons
		G4StepPoint*      startpoint = step->GetPreStepPoint();
		G4TouchableHandle geometry   = startpoint->GetTouchableHandle();
		G4String          volumeName = geometry->GetVolume()->GetName();
		static bool       first      = true;
		static int        up         = 0;
		static int        down       = 0;
		static G4double     quantum_eff = 1.;

		bool record_hit = (G4UniformRand()<quantum_eff);

		if (volumePositions.count(volumeName) > 0){ 
			if (record_hit && (volumeName=="up" || volumeName =="down")){
				if (first){ G4cout << "#DET\tdetector\ttime[ns]\tOrigin Position (x,y,z)\tTrackLen" << G4endl; first = false;}
				if (volumeName == "up"){up++;}
				if (volumeName == "down"){down++;}
				G4double time = startpoint->GetGlobalTime()/ns;
				//G4ThreeVector dir = -startpoint->GetMomentumDirection();
				G4ThreeVector createdPos = step->GetTrack()->GetVertexPosition();
				G4double tracklen = createdPos.mag();
				//dir /= dir.mag();
				G4cout << "DET" << "\t" << volumeName << "\t" << time  << "\t" << createdPos << "\t" << tracklen << G4endl;
			}
			if (volumeName=="up" || volumeName == "down" || volumeName == "housing"){
				step->GetTrack()->SetTrackStatus(fStopAndKill);
			}
		}

		if (step->GetTrack()->GetParentID()==0 && step->GetTrack()->GetCurrentStepNumber()==1 && !(first)){
			if (up>0 && down>0) {
				int coincCount = up<down ? up : down;
				G4cout << "Coincidences:\t" << coincCount << G4endl;
			}
			G4cout << "\nNew Track!" << G4endl;
			first = true;
			up = 0; down = 0;

		}
		
	}
};
		
	
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
		u->ApplyCommand("/tracking/verbose 1");
		u->ApplyCommand(G4String("/run/beamOn ") +G4String(to_string((int)var["n"])));
		ui->SessionStart();
		delete ui;
	}
	else if (var["Verbose"] != 0.0){
		u->ApplyCommand("/tracking/verbose 1");
		u->ApplyCommand(G4String("/run/beamOn ") +G4String(to_string((int)var["n"])));
	}
	else{
		u->ApplyCommand(G4String("/run/beamOn ") +G4String(to_string((int)var["n"])));
	}
	delete visManager;
	delete runManager;
}
