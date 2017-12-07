
/*______________________________________________________________________________
|                                                                              |
|o------------------------------------>X<-------------------------------------o|
|                                                                              |
| ___________________________________________________________________________  |
| |   This program was written by:           4/28/2014                       | |
| |    Samuel David Lotz                                                     | |
| |     Title:                                                               | |
| |	 Protein_Folder                                                          | |
| |       Protein Folder is a project written under the guidance of
| |   Dr. Sam Thangiah at Slippery Rock University for a Machine learning 
| |   independent instruction course.  Code design and pseudo-
| |   code taken from Alena Shmygelska's dissertation "Novel Heuristic
| |   Search Methods for Protein Folding and Identification of Folding 
| |   Pathways" Sept 2006 University of British Columbia. The goal of this
| |   project was to implement machine learning strategies to the 
| |   biological problem of protein folding and to compare them to non-
| |   Learning stategies.
| | 
| |__________________________________________________________________________| |
|______________________________________________________________________________|*/
/* -------   Class-------------------------------------------------------

	 
	
	ATTRIBUTES:
	
	
	CONSTRUCTORS:

	
	METHODS:
	GETTERS
	SETTERS
	
	
	
	
	Version 1.0 beta 
		
  ------------------------------------------------------------------------------
 */


package protein_folder;

import java.awt.Color;
import java.util.ArrayList;

public class Param{

//==================Standard Constants======================================================
	final public static double ZERO = 0;
	final public static double UNITY = 1.0;
	final public static double THOUSAND = 1000.0;
	final public static int LATTICE_UNIT = 1;  //the standard length of a residue bond length
//==================================================================================	
	
//==================Enumerated Types======================================================
	
	public static enum Direction {S, L, R, E} //straight ahead, left, right, and end of chain
	public static enum ResidueType {H , P, NULL} //H = hydrophobic, P = polar
	public static enum TerminalType {NOT, C, N} //for not a terminal residue c-terminus, and N-terminus
	public static enum UnitVector {i, j, I, J, ZERO} //the directions in which the residues can be oriented towards
		//lower case is positive direction capital is negative direction
		//named unit vectors because we may use other lattices in the future
		//i is +x, I is -x, j is +y, J is -y
	public static enum Options {RANDOM, PHEROMONE_TRAILS, LOCAL_SEARCH, PROTEIN_BUILD}  //allows for the activation of different 
			//machine learning algorithm types for a colony
//==================================================================================
	
//==================Output Constants======================================================

	
	public static int FONT_SIZE = 20;
	public static Color FONT_COLOR = Color.BLACK;
	public static int MARGINS = 20;
	public static int LATTICE_SIZE = 100;
	public static int FOOTER_SIZE = 200;
	public static int RESIDUE_CIRCLE_DIAMETER = 50;
	public static int RESIDUE_CIRCLE_STROKE_WIDTH = 3;
	public static Color POLAR_COLOR = Color.WHITE;
	public static Color HYDROPHOBIC_COLOR = Color.BLUE;
	public static int CONNECTOR_WIDTH = 15;
	public static Color CONNECTOR_COLOR = Color.BLACK;
	public static Color BACKGROUND_COLOR = Color.WHITE;
	public static Color STROKE_COLOR = Color.BLACK;
	

//==================================================================================

//==================Enum Type Public Keywords============================================
	//use these for the enum types instead of calling enumType.Type use Param.keyword
	final public static Direction[] allDirections2D = {Direction.S, Direction.L, Direction.R};
	//final public static Direction[] allDirections3D = {Direction.U, Direction.D, Direction.S, Direction.L, Direction.R}; not functional yet in 3D	
	final public static Direction STRAIGHT = Direction.S;
	final public static Direction LEFT = Direction.L;
	final public static Direction RIGHT = Direction.R;
	final public static Direction END = Direction.E;
	final public static ResidueType POLAR = ResidueType.P;
	final public static ResidueType HYDROPHOBIC = ResidueType.H;
	final public static ResidueType NULL = ResidueType.NULL;
	final public static TerminalType C_TERMINUS = TerminalType.C;
	final public static TerminalType N_TERMINUS = TerminalType.N;
	final public static TerminalType NOT_TERMINUS = TerminalType.NOT;
	final public static UnitVector i_VECTOR = UnitVector.i;
	final public static UnitVector I_VECTOR = UnitVector.I;
	final public static UnitVector j_VECTOR = UnitVector.j;
	final public static UnitVector J_VECTOR = UnitVector.J;
	final public static UnitVector ZERO_VECTOR = UnitVector.ZERO;
	final public static Options RAND_OPT = Options.RANDOM;
	final public static Options PHER_OPT = Options.PHEROMONE_TRAILS;
	final public static Options LOC_SEARCH_OPT = Options.LOCAL_SEARCH;
	
//==================================================================================
	
//==================Default Parameters======================================================
	//From Shmygelska dissertation
	final public static double DEF_ALPHA = 1.0;
	final public static double DEF_BETA = 2.0;
	final public static double DEF_RHO = 0.8;
	final public static double DEF_THETA = 0.05;
	final public static double DEF_INV_TEMP = 26;
	final public static double INIT_PHEROMONE_2D = 0.33333333333;
	
	final public static double DEF_SELECTION_THRESHOLD = 0.1;    //should be 1% (.01) according to Shmygelska
	final public static double MIN_SEL_THRESH = 0.01; //1% FOR 100 ANTS IS 1 ANT PER WAVE ANY LESS NO UPDATE WOULD OCCUR
	final public static int DEF_LOCAL_SEARCH_MAX_ITER_N_LESS_50 = 1000;
	final public static int DEF_LOCAL_SEARCH_MAX_ITER_N_GREATER_50 = 10000;
	final public static double DEF_KEEP_CONF_LOCAL_SEARCH = 0.5;
	final public static double DEF_IMPROVING_ANTS = 0.5;  //the proportion of ants in the colony which are not allowed to local search
	
	final public static int DIRECTIONS_2D = 3;
	
//==================================================================================
	
//==================Sequence Library======================================================
	final public static int NUM_SEQ_LIB = 12;
	final public static int[] SEQ_LIBRARY = {0,1,2,3,4,5,6,7,8,9,10,11};
//from Shmygelska dissertation
	//default sequence
	final private static ResidueType[] defaultSeq = {POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR};
	final private static int defSeq_Length = 50;
	final private static double defSeq_OptE = -28;
	final private static double seq0_cutoff = -28;
	//2D standard sequences
		//sequence 1
	final private static ResidueType[] seq1 = {HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC}; 
	final private static int seq1_Length = 20;
	final private static double seq1_OptE = -9;
	final private static double seq1_cutoff = -9;
		//seuence 2
	final private static ResidueType[] seq2 = {HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC};
	final private static int seq2_Length = 24;
	final private static double seq2_OptE = -9;
	final private static double seq2_cutoff = -9;
		
		//seuence 3
	final private static ResidueType[] seq3 = {POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC};
	final private static int seq3_Length = 25;
	final private static double seq3_OptE = -8;
	final private static double seq3_cutoff = -8;
		//seuence 4
	final private static ResidueType[] seq4 = {POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR};
	final private static int seq4_Length = 36;
	final private static double seq4_OptE = -14;
	final private static double seq4_cutoff = -14;
		//seuence 5
	final private static ResidueType[] seq5 = {POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC};
	final private static int seq5_Length = 48;
	final private static double seq5_OptE = -23;
	final private static double seq5_cutoff = -23;
		//seuence 6
	final private static ResidueType[] seq6 = {HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC};
	final private static int seq6_Length = 50;
	final private static double seq6_OptE = -21;
	final private static double seq6_cutoff = -21;
		//seuence 7
	final private static ResidueType[] seq7 = {POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR};
	final private static int seq7_Length = 60;
	final private static double seq7_OptE = -36;
	final private static double seq7_cutoff = -21;
		//seuence 8
	final private static ResidueType[] seq8 = {HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,POLAR,HYDROPHOBIC,POLAR,HYDROPHOBIC,POLAR,POLAR,HYDROPHOBIC,HYDROPHOBIC,POLAR,POLAR,HYDROPHOBIC,HYDROPHOBIC,POLAR,POLAR,HYDROPHOBIC,POLAR,POLAR,HYDROPHOBIC,HYDROPHOBIC,POLAR,POLAR,HYDROPHOBIC,HYDROPHOBIC,POLAR,POLAR,HYDROPHOBIC,POLAR,POLAR,HYDROPHOBIC,HYDROPHOBIC,POLAR,POLAR,HYDROPHOBIC,HYDROPHOBIC,POLAR,POLAR,HYDROPHOBIC,POLAR,HYDROPHOBIC,POLAR,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC,HYDROPHOBIC};
	final private static int seq8_Length = 64;
	final private static double seq8_OptE = -42;
	final private static double seq8_cutoff = -42;
		//seuence 9
	final private static ResidueType[] seq9 = {HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC};
	final private static int seq9_Length = 85;
	final private static double seq9_OptE = -53;
	final private static double seq9_cutoff = -53;
		//seuence 10
	final private static ResidueType[] seq10 = {POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC};
	final private static int seq10_Length = 100;
	final private static double seq10_OptE = -50;
	final private static double seq10_cutoff = -50;
		//seuence 11
	final private static ResidueType[] seq11 = {POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC, HYDROPHOBIC, POLAR, POLAR, POLAR, POLAR, POLAR, POLAR, HYDROPHOBIC, POLAR, HYDROPHOBIC, HYDROPHOBIC};
	final private static int seq11_Length = 100;
	final private static double seq11_OptE = -48;
	final private static double seq11_cutoff = -48;
//==================================================================================
	
//==================ATTRIBUTES======================================================
	//Experiment Parameter attributes
	private String title = "TITLE NOT INIT";

	private ArrayList<Options> option = new ArrayList<Options>(); //what options are activated for this colony
	private int numDirections;
	private ArrayList<Direction> allDirections = new ArrayList<Direction>();
	private double invTemp; //the inverse temperature used in calculating the eta values for new folds
	private double alpha; // a parameter used to weight the pheromone weight in new fold decisions
	private double beta; // a parameter used to weight the eta weight in new fold decisions
	private double rho; //the evaporation rate of pheromone
	private double theta; //the minimal pheromone value for renormalization in the MAX-MIN ant system
	private double initPheromone; //the initial amount of pheromone put down on the trails
	private double selectionThreshold; //the percentage of best conformations that will update the pheromone matrix
	private int localSearchMaxIter; //the number of attempts the local search will make before it stops if no better solution is found
	private double probKeepConfLocalSearch;  //the probability the previous conformation in a long range local search is kept, when feasible
	private double proportionImprovingAnts;  ////the proportion of ants in the colony which are not allowed to local search
	

	//sequence specific parameters
	private ArrayList<ResidueType> residueSeq = new ArrayList<ResidueType>(); //array of the residue types to fold, input target sequence
	private String sequenceName;
	private int proteinLength; //the number of residues in the protein being folded
	private double optimalEnergy; //the theoretical optimal energy that will be compared to in the pheromone update
	private double solutionCutoff;
	
//==================================================================================
		
//////__________Constructors_________________________________________________________

//----No-Arg Constructor-----------------------------------------------------------
	//makes default parameters for a 2D protein
	public Param(String inTitle){
		/*
		//set the parameters
				title = inTitle;
				sequenceName = "Default_2D";
				option.add( RAND_OPT);
				numDirections = DIRECTIONS_2D;
				for(int i=0; i < numDirections; i++){
					allDirections.add(allDirections2D[i]);
				}
				invTemp = 1;
				alpha = DEF_ALPHA;
				beta = DEF_BETA;
				rho = DEF_RHO;
				theta = DEF_THETA;
				initPheromone = INIT_PHEROMONE_2D;
				selectionThreshold = DEF_SELECTION_THRESHOLD;
				proteinLength = defSeq_Length;
				setOptimalEnergy(0);
				for(int i=0; i<proteinLength; i++){
					residueSeq.add(seq1[i]);
				}
				*/
	}//end default no-arg constructor
//----End No-Arg Constructor-------------------------------------------------------
	
//----Default Constructor-----------------------------------------------------------
	// sets default parameters for a sequence from the library
	public Param( int seqNum, String inTitle){
		
		//set the parameters
		title = inTitle;
		sequenceName = "Sequence_" + seqNum + "_2D";
		option.add(RAND_OPT);
		numDirections = DIRECTIONS_2D;
		for(int i=0; i < numDirections; i++){
			allDirections.add(allDirections2D[i]);
		}
		invTemp = DEF_INV_TEMP;
		alpha = DEF_ALPHA;
		beta = DEF_BETA;
		rho = DEF_RHO;
		theta = DEF_THETA;
		initPheromone = INIT_PHEROMONE_2D;
		selectionThreshold = DEF_SELECTION_THRESHOLD;
		
		//now set the sequence specific parameters
		this.setSeqInfo(seqNum);
		
		//set local search parameters
		probKeepConfLocalSearch = DEF_KEEP_CONF_LOCAL_SEARCH;
		proportionImprovingAnts = DEF_IMPROVING_ANTS;
		
		//set the localSearchMaxIter for proteins of different sizes
		if(this.proteinLength <= 50)
			localSearchMaxIter = DEF_LOCAL_SEARCH_MAX_ITER_N_LESS_50;
		else if(this.proteinLength > 50)
			localSearchMaxIter = DEF_LOCAL_SEARCH_MAX_ITER_N_GREATER_50;
		
	}//end default arg constructor
//----End Default Constructor-------------------------------------------------------
	
//___________________________________________________________________________________	
	
	
	
//----setSeqInfo-----------------------------------------------------------
	
	//when this is called it sets the current fields for the Param object to the specified sequence
	private void setSeqInfo(int inSeqNum){
			switch(inSeqNum){

			case 0: proteinLength = defSeq_Length;
					setOptimalEnergy(defSeq_OptE);
					solutionCutoff = seq0_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(defaultSeq[i]);
					}
				break;
				
			case 1: proteinLength = seq1_Length;
					setOptimalEnergy(seq1_OptE);
					solutionCutoff = seq1_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(seq1[i]);
					}
				break;
				
			case 2: proteinLength = seq2_Length;
					setOptimalEnergy(seq2_OptE);
					solutionCutoff = seq2_cutoff;
						for(int i=0; i<proteinLength; i++){
								residueSeq.add(seq2[i]);
						}
				break;
			case 3: proteinLength = seq3_Length;
					setOptimalEnergy(seq3_OptE);
					solutionCutoff = seq3_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(seq3[i]);
					}
				break;
			case 4: proteinLength = seq4_Length;
					setOptimalEnergy(seq4_OptE);
					solutionCutoff = seq4_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(seq4[i]);
					}
				break;
			case 5: proteinLength = seq5_Length;
					setOptimalEnergy(seq5_OptE);
					solutionCutoff = seq5_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(seq5[i]);
					}
				break;
			case 6: proteinLength = seq6_Length;
					setOptimalEnergy(seq6_OptE);
					solutionCutoff = seq6_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(seq6[i]);
					}
				break;
			case 7: proteinLength = seq7_Length;
					setOptimalEnergy(seq7_OptE);
					solutionCutoff = seq7_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(seq7[i]);
					}
				break;
			case 8: proteinLength = seq8_Length;
					setOptimalEnergy(seq8_OptE);
					solutionCutoff = seq8_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(seq8[i]);
					}
				break;
			case 9: proteinLength = seq9_Length;
					setOptimalEnergy(seq9_OptE);
					solutionCutoff = seq9_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(seq9[i]);
					}
				break;
			case 10: proteinLength = seq10_Length;
					setOptimalEnergy(seq10_OptE);
					solutionCutoff = seq10_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(seq10[i]);
					}
				break;
			case 11: proteinLength = seq11_Length;
					setOptimalEnergy(seq11_OptE);
					solutionCutoff = seq10_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(seq11[i]);
					}
				break;
			default: proteinLength = defSeq_Length;
					setOptimalEnergy(defSeq_OptE);
					solutionCutoff = seq0_cutoff;
					for(int i=0; i<proteinLength; i++){
						residueSeq.add(defaultSeq[i]);
					}
				break;
		
		}//end switch
			
	}//end the setSeqInfo
	
//----End setSeqInfo-------------------------------------------------------
	
	
	
	
//////++++++++GET String Methods+++++++++++++++++++++++++++++++++++++++++++++++++++++
	public String getStringDirection(Direction in){
		String out = null;
		switch(in){
			case S: out = "S";
				break;
			case L: out = "L";
				break;
			case R: out = "R";
				break;
			default: out = "E";
				break;
		}//end switch
		return out;
	}//end getCharDirection
	
	public String getStringResType(ResidueType in){
		String out = null;
		switch(in){
			case H: out = "H";
				break;
			case P: out = "P";
				break;
			default: out = "NULL";
				break;
		}//end switch
		return out;
	}//end getCharResType
	
	public String getStringUnitVector(UnitVector in){
		String out = null;
		switch(in){
			case i: out = "i";
				break;
			case I: out = "I";
				break;
			case j: out = "j";
				break;
			case J: out = "J";
				break;
			default:out = "NULL";
				break;
		}//end switch
		return out;
	}//end getStringUnitVector
	
	public String getStringTerminus(TerminalType in){
		String out = null;
		switch(in){
			case N: out = "N";
				break;
			case C: out = "C";
				break;
			case NOT: out = "NOT";
				break;
			default: out = "NULL";
						break;
		}//end switch
		return out;
	}//end getCharTerminus
	
//////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	
//////++++++++GETTERS & SETTERS+++++++++++++++++++++++++++++++++++++++++++++++++++++

	public double getSolutionCutoff() {
		return solutionCutoff;
	}

	public void setSolutionCutoff(double solutionCutoff) {
		this.solutionCutoff = solutionCutoff;
	}

	public double getProportionImprovingAnts() {
		return proportionImprovingAnts;
	}
	public void setProportionImprovingAnts(double proportionLocSearchAnts) {
		this.proportionImprovingAnts = proportionLocSearchAnts;
	}
	public double getProbKeepConfLocalSearch() {
		return probKeepConfLocalSearch;
	}
	public void setProbKeepConfLocalSearch(double probKeepConfLocalSearch) {
		this.probKeepConfLocalSearch = probKeepConfLocalSearch;
	}
	public int getLocalSearchMaxIter() {
		return localSearchMaxIter;
	}
	public void setLocalSearchMaxIter(int localSearchMaxIter) {
		this.localSearchMaxIter = localSearchMaxIter;
	}
	public String getSequenceName() {
		return sequenceName;
	}
	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}
	public String getTitle() {
		return title;
	}
	public void setTitle(String title) {
		this.title = title;
	}
	public ArrayList<Options> getOption() {
		return option;
	}
	public void setOption(ArrayList<Options> option) {
		this.option = option;
	}
	public void setAllDirections(ArrayList<Direction> allDirections) {
		this.allDirections = allDirections;
	}
	public ArrayList<Direction> getAllDirections() {
		return allDirections;
	}
	public int getNumDirections() {
		return numDirections;
	}
	public void setNumDirections(int numDirections) {
		this.numDirections = numDirections;
	}
	public double getInvTemp() {
		return invTemp;
	}
	public void setInvTemp(double invTemp) {
		this.invTemp = invTemp;
	}
	public double getAlpha() {
		return alpha;
	}
	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}
	public double getBeta() {
		return beta;
	}
	public void setBeta(double beta) {
		this.beta = beta;
	}
	public double getRho() {
		return rho;
	}
	public void setRho(double rho) {
		this.rho = rho;
	}
	public double getTheta() {
		return theta;
	}
	public void setTheta(double theta) {
		this.theta = theta;
	}
	public ArrayList<ResidueType> getResidueSeq() {
		return residueSeq;
	}
	public void setResidueSeq(ArrayList<ResidueType> residueSeq) {
		this.residueSeq = residueSeq;
	}
	public double getProteinLength() {
		return proteinLength;
	}
	public void setProteinLength(int proteinLength) {
		this.proteinLength = proteinLength;
	}
	public double getSelectionThreshold() {
		return selectionThreshold;
	}
	public void setSelectionThreshold(double selectionThreshold) {
		this.selectionThreshold = selectionThreshold;
	}
	public double getOptimalEnergy() {
		return optimalEnergy;
	}
	public void setOptimalEnergy(double optimalEnergy) {
		this.optimalEnergy = optimalEnergy;
	}
	public double getInitPheromone() {
		return initPheromone;
	}
	public void setInitPheromone(double initPheromone) {
		this.initPheromone = initPheromone;
	}

//////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



	
}//end param class