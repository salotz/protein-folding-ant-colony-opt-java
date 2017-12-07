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

import static java.nio.file.StandardOpenOption.APPEND;
import static java.nio.file.StandardOpenOption.CREATE;
import static java.nio.file.StandardOpenOption.WRITE;

import java.awt.image.BufferedImage;
import java.io.*;
//data structures
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Random;
import java.lang.StringBuilder;
import java.lang.Math;
import java.nio.charset.Charset;
import java.nio.file.*;

import javax.imageio.ImageIO;

import protein_folder.Param.Direction;
import protein_folder.Param.ResidueType;
import protein_folder.Param.UnitVector;
import protein_folder.Param.Options;

public class Driver {
	
	//TO THE USER
	//SET THE EXPERIMENT SPECIFICS
	//default published way to do this is
	//100 ANTS
	//for seqLength <= 50 500 waves,  50 < seqLength <= 64 100 waves,  seqLength < 64 20 waves
	
	static String EXPERIMENT_BIN = "C:\\Users\\Sam Lotz\\Desktop\\Protein_Folder Data\\" ;  //the folder where the experiment folder will go
	static int NUM_ANTS = 100;
	static int NUM_WAVES = 500;
	//static int NUM_WAVES = 100;
	//static int NUM_WAVES = 20;
	static int[] SEQUENCES = {1};
	
	static int buildSeqNum = 10; //was using this to get a few Protein_Build objects saved and printed.  Doesn't really work not important.
	
	public static void main(String[] args)
	{
		double startTime = System.nanoTime();
		
		//create a directory for all files to go in for this experiment
		Path expDir;
		//set the experiment directory path
		
		try {
			expDir = Paths.get( initExpDir(setExpTitleDateTime()).toString() );
			//if the name already exists exit program and try again
		} catch (FileAlreadyExistsException e1) {
			
			//SYSTEM OUTPUT
			System.out.println("Aborting Experiment");
			return ;
		}//end catch
		
		
		//make parameter
		Param param = new Param(1, "seq1_def_param");
		
		//make the protein with first residue
		Protein prot = new Protein(param.getResidueSeq().get(0));
		
		//create a ProteinBuild and add the first protein to it
		ProteinBuild protBuild = new ProteinBuild();
		protBuild.getBuildSteps().add(new Protein(prot));
		
		//build a straight chain version of the protein and add each step to the build
		for(int i=1; i < param.getProteinLength(); i++){
			prot.addResidue(Param.STRAIGHT, param.getResidueSeq().get(i));
			protBuild.getBuildSteps().add(new Protein(prot));
		}
		
		
		
		
		//make images from the build
		PNG_Footer footer = new PNG_Footer();
		
		for(int i=0; i < protBuild.getBuildSteps().size(); i++){
			
			//set footer specifics
			footer.setSolutionStr("Step :  " + i) ;
			footer.setEnergyStr("Energy:  " + protBuild.getBuildSteps().get(i).getEnergy());
			footer.setSequenceInfoStr("Sequence: "+param.getSequenceName() + "  Length: "+ param.getProteinLength()+ "  Opt E:  "+ param.getOptimalEnergy());
			footer.setColInfoStr("Colony: " );
			
			try {
				protBuild.getBuildStepsImages().add( Protein.makeImage( protBuild.getBuildSteps().get(i), footer ) );
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		
		
		}//end make images for build
		
		
		try {
			Path buildDir;
			//make a title for this build's directory
			String buildTitle = "Build";
			buildDir = initBuildDir(expDir, buildTitle);
			
			ListIterator<BufferedImage> imgsIter = protBuild.getBuildStepsImages().listIterator();
			int imgCount = 1;
			while(imgsIter.hasNext()){
				
				ImageIO.write( imgsIter.next(),  "PNG", new File(buildDir.toString() + "\\step_" + imgCount + ".PNG") );
				imgCount++;
			}
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
		
		
		
		
		
		
		
		
		/*
		
		//DEBUGGING experimental set-up
		ArrayList<Colony> expCols = new ArrayList<Colony>();
		Param param1 = new Param(1, "defParam_noPher_seq1");
		expCols.add( new Colony(param1, 1, 1, "debug_colony_waves1_ants1") );
		*/
		
		
		
		/*
		//Experiment set-up
		//ArrayList<Param> expParamsNoPher = new ArrayList<Param>();
		ArrayList<Param> expParamsPher = new ArrayList<Param>();
		//ArrayList<Param> expParamsPherLocal = new ArrayList<Param>();
		
		//define parameter sets to be used
		for(int i=0; i< SEQUENCES.length; i++){
			
			//No Pheromone Params
			
		//	expParamsNoPher.add( new Param( SEQUENCES[i], "defParam_noPher_seq"+SEQUENCES[i]) );
			
			//SYSTEM OUTPUT
			System.out.println("No Pheromone Param Created");
			
			expParamsPher.add(  new Param( SEQUENCES[i], "defParam_Pher_seq"+SEQUENCES[i]) );
			expParamsPher.get(i).getOption().add(Param.PHER_OPT);
			
			//SYSTEM OUTPUT
				System.out.println("Pheromone Param Created");
				
				
			//expParamsPherLocal.add(  new Param( SEQUENCES[i], "defParam_PherLocal_seq"+SEQUENCES[i]) );
			//expParamsPherLocal.get(i).getOption().add(Param.PHER_OPT);
			//expParamsPherLocal.get(i).getOption().add(Param.LOC_SEARCH_OPT);
			//SYSTEM OUTPUT
				System.out.println("Pheromone,Local Search Param Created");
				
		}//end create Params for each sequence
				
		
		ArrayList<Colony> expCols = new ArrayList<Colony>();
		//colony constructor (param, numAnts, waves)
		
		
		//make colonies without pheromone
		for(int i=0; i < expParamsNoPher.size(); i++){
			
			//SYSTEM OUTPUT
			System.out.println("No Pheromone Colony "+ i + ":");
			expCols.add( new Colony(expParamsNoPher.get(i), NUM_ANTS, NUM_WAVES, NUM_ANTS + "ant_by_" + NUM_WAVES + "_waves_" +  expParamsNoPher.get(i).getTitle() + "_" + i) );
			
			//SYSTEM OUTPUT
			System.out.println("No Pheromone Colony "+ i + " created");
		}//end making no pheromone colonies
		
		
		//make colonies with pheromone
		for(int i=0; i < expParamsPher.size(); i++){
			//SYSTEM OUTPUT
			System.out.println("Pheromone Colony "+ i + ":");
			
			expCols.add( new Colony(expParamsPher.get(i), NUM_ANTS, NUM_WAVES, NUM_ANTS + "ant_by_" + NUM_WAVES + "_waves_" + expParamsPher.get(i).getTitle() + "_" + i) );
			
			//SYSTEM OUTPUT
			System.out.println("Pheromone Colony "+ i + " created");
		}//end making pheromone trail sequences
		
		
		//make colonies pheromone and Local searching
		for(int i=0; i < expParamsPherLocal.size(); i++){
			//SYSTEM OUTPUT
			System.out.println("Pheromone,Local Search Colony "+ i + ":");
			
			expCols.add( new Colony(expParamsPherLocal.get(i), NUM_ANTS, NUM_WAVES, NUM_ANTS + "ant_by_" + NUM_WAVES + "_waves_" + expParamsPherLocal.get(i).getTitle() + "_" + i ) );
			
			//SYSTEM OUTPUT
			System.out.println("Pheromone, Local Search Colony "+ i + " created");
		}//end making pheromone trail sequences
		
		
		*/
		
		/*
	//FILE OUTPUT
		//each colonies output
		ListIterator<Colony> colIter = expCols.listIterator();
		//iterate through the list of colonies and write a file for each
		while(colIter.hasNext()){
			writeColonyDataFile(colIter.next(), expDir);
			
			//SYSTEM OUTPUT
			System.out.println("Colony files created");
		}
		
		
		//output for comparing wave averages between the colonies
		writeTotalResults(expDir, expCols);
		
		*/
	/*//make build sequence images and directories
		makeBuildSeqs(expCols, expDir);*/
		
		
		
	double endTime = System.nanoTime();
	double totalTime = endTime-startTime;
	
	//SYSTEM OUTPUT
	System.out.println("Total Run Time: " + totalTime);
		
	}//end main
	
	//----localSearch-----------------------------------------------------------

		//takes as an argument the proteins current conformation
		//returns either the same protein conformation or one with an improved energy conformation
		//will attempt to make better conformations either until it finds a better one
		//of the localSearchMaxIter number has been reached and none have been found
		//the probability that a residue changes its relative next direction is given in the colParam
		static public Protein localSearch(Protein currConf, Param inParam){
			
			Protein outProt = new Protein(); //the protein the function will return
			Random rand = new Random();
			int mutResNum = 0;  // the residue index that the local search long range movement starts at
			boolean feasible = false; //flag raised when the current residues direction is feasible
			Direction tempDir = Param.END; 
			boolean abort = false; //raised when the current local search is problematic and should be aborted
			boolean solution = false; //raised when a conformation is found that is better than the current one
			
			int i=0;
			//create new conformations until a better one is found or the limit for searches is up
			while(i < inParam.getLocalSearchMaxIter()  &&   !solution){
				//make a copy of this protein to manipulate
				Protein newConf = new Protein(currConf);
				
				//randomly choose a residue to change
				//dont choose the first residue because we can not calculate last unit vector
				mutResNum = 1 + rand.nextInt(currConf.getLength()-1);
				
				//remove the residues after this one
				for(int j=1; j < (currConf.getLength() - mutResNum); j++){
					newConf.removeResidue();
				}//end removing residues after the ith chosen one for the new conformation
				newConf.updateMap();
				
				
	//DEBUG*****************************************
				/*System.out.println("The current conformation");
				currConf.printCoordMap();
				System.out.println("The new conformation");
				newConf.printCoordMap();*/
				
				
				//now go through the rest of the residues and change their directions is either the direction is not feasible anymore
				//of  with probability = ( 1 - probability of keeping current change)
				for(int j=0; j < (currConf.getLength() - mutResNum); j++){
					
					//for this residue find the possible directions
					ArrayList<Direction> feasibleNewDirs = new ArrayList<Direction>();
					feasibleNewDirs.addAll( findPossDirections(newConf, inParam) );
					
					//check to see if the currConf previous position is infeasible
					feasible = false; //if a match is found in possDirs then it is
					for(int k=0; k < feasibleNewDirs.size(); k++){
						if(currConf.getSequence().get(mutResNum + j).getNextResDirection() == feasibleNewDirs.get(k)){
							feasible = true;
						}
					}//end check for feasibility
					
					//if probability wills it or the previous direction is infeasible make new direction
					if( rand.nextDouble() <= inParam.getProbKeepConfLocalSearch() ||  !feasible) {
						
						//set the proposed fold to tempDirection
						tempDir = addFoldLocalSearch(currConf.getSequence().get(mutResNum + j).getType(), newConf, inParam);
						
						//ensure that there is a feasible direction
						if(tempDir == Param.END){
							//need to start the local search over
							i--;
							abort = true;
						}else{
							//there is a possible direction for the new residue so add it
							newConf.addResidue(tempDir, currConf.getSequence().get(mutResNum+j).getType() );
							
						}//end no possDirectinos guard
						
					}//end make new direction guard
					
					//if the current local search is problematic break out of the j rebuild loop
					//reseed the local search
					if(abort)
						break;
					
				}//end change residues for this local search
				
				//check to see if the resulting long range move is an improvement
				if(newConf.getEnergy() < currConf.getEnergy()){
					//if it is better stop looking and set the output protein to the newConf
					solution= true;
					outProt = newConf;
				}
					
				
				
			}//end make new possible conformations
			
			//if a better solution was found return it
			if(solution){
				return outProt;
			}
			else{ //if not give back the old one
				return currConf;
			}
			
			
		}//end localSearch
		
	//----End localSearch-------------------------------------------------------
		
		//----addFoldLocalSearch-----------------------------------------------------------	
		//adds a fold to the protein by finding the options, making a choice, and then returning the direction
		// this is a large function with multiple phases which I have included all in one for now because the sub-functions
		//require a lot of parameters and it might be easier to house them under one function
		//Phase 1: find all the spots that the new residue could possibly be put
		//Phase 2: if no spots do attrition
		//Phase 3: Make a probabilistic choice based on the pheromone and thermo favorability
		//returns the direction of the next residue
		static public Direction addFoldLocalSearch(ResidueType inType, Protein inProt, Param inParam){
			Random rand = new Random();
			//the direction we will return
			Direction newDirection = Param.END;
			double choiceNum = rand.nextDouble(); //used to pick a direction (0.0 <= choiceNum < 1.0)
			double cumProb = 0; //used for calculating cumulative probabilities
			
			//make a list of the possible directions and fill it
			ArrayList<Direction> possDirections = new ArrayList<Direction>();
			possDirections.addAll(findPossDirections(inProt, inParam));
			
			//the cumulative distribution boundaries, first element is 0
			Double[] probs = new Double[possDirections.size() + 1];
			probs[0] = cumProb; //because it is 0 right now
			
			
			//if there are no possible directions we must start the local search over
			//so return the null direction to signal to the local search to start over
			if( possDirections.isEmpty()){
				
	//DEBUG***********************************************			
				//System.out.println("no possible directions local search begins again");
				
				
				
				return Param.END;
			}
			else{// we choose which direction to return based on the thermodynamic model
			
				//first create the cumulative probability distribution
				cumProb=0;
				for(int i=0; i < possDirections.size(); i++){
					probs[i+1] = (cumProb + calcFoldProbLocalSearch(inProt, inType, possDirections.get(i), inParam));  //DEBUG take out the first parameter
					//increase the cumulative probability
					cumProb = probs[i+1]; //because the first element is 0
	/*				
			//DEBUG
					System.out.println("Cumulative probability " + i + " : " + cumProb);
					*/
				}//end for 
				
				//Then pick, use the cumulative distribution to find which direction is the choice
				int i=0;
				//while the newDirection has not been set to a direction or the end of the possible directions
				while( newDirection == Param.END && i< possDirections.size()){
					
					//check to see which part of the inequality the chosen valuse lies
					if(probs[i] <= choiceNum && choiceNum < probs[i+1] ){
						//then set the newDirection to the this direction
						newDirection = possDirections.get(i);
						
					}//end the check for a choice and newDirection set
					i++;
				}//end while search through possible directions 
			
			return newDirection;
			}//end if else
		}//end addFoldLocalSearch method
	//----End addFoldLocalSearch-------------------------------------------------------
		
		//----calcFoldProbLocalSearch-----------------------------------------------------------
		//the calculate fold probability method
		public static double calcFoldProbLocalSearch(Protein inProt, ResidueType inType, Direction foldDirection, Param inParam){
			double outProb = 0.0;
			double eta = 0;
			int newHH = 0; //the new H-H bonds created by the input fold
			double weightNumerator = 0.0; //the weight of the fold of interest i.e. the numerator of the probability function
			double weightTotal = 0.0;  //the denominator in the probability function
		
			//first find the numerator or the value of the fold in question
			
			//protein that is constructed to calculate the energy from the input fold
			Protein testProtNum = new Protein(inProt);
			//array for possible folds
			ArrayList<Direction> possFolds = new ArrayList<Direction>();
			
			//calculate the number of new H-H bonds by making the test protein and calculating the energy difference
			testProtNum.addResidue(foldDirection, inType);
			// the getEnergy automatically updates the energy -- testProt.updateEnergy();
			newHH = testProtNum.getEnergy() - inProt.getEnergy();
			//calculate the eta value
			eta = Math.exp(-inParam.getInvTemp() * newHH);
			//find the tau value from the pheroMatrix

			//save to the numerator value
			weightNumerator = eta;
			
		//DEBUG*********************************************
			//System.out.println("The numerator : " + weightNumerator);
			
			
			//NEXT SUM OF ALL POSSIBLE FOLDS
			//calculate the sum of all the possible folds
			possFolds = findPossDirections(inProt, inParam);
			//calculate the weights and add them up
			for(int i=0; i<possFolds.size(); i++){
								
				Protein testProtTemp = new Protein(inProt);
				//calculate the number of new H-H bonds by making the test protein and calculating the energy difference
				testProtTemp.addResidue(possFolds.get(i), inType);
				// the getEnergy automatically updates the energy -- testProt.updateEnergy();
				newHH = testProtTemp.getEnergy() - inProt.getEnergy();
				//calculate the eta value
				eta = Math.exp(-inParam.getInvTemp() * newHH);
				//find the tau value
				//STUBBED OUT
				
				//add it to the weight total
				weightTotal = weightTotal +eta;
				
		//DEBUG**************************************************
				//System.out.println("Denominator term "+ i + weightTotal);
			
			}//end for totalling weights

		//DEBUG*******************************************************
			//System.out.println("The denominator : " + weightTotal);
			
			//calculate the probability
			outProb = weightNumerator/weightTotal ;
			
		//DEBUG************************************************************
			//System.out.println("The probability of fold: " + outProb);
			
			
			return outProb;
		}//end calcFoldProbLocalSearch
	//----End calcFoldProbLocalSearch-------------------------------------------------------
		
		//----getLastUnitVector-----------------------------------------------------------
		//get the unit vector direction to the last residue
		public static UnitVector getLastVectorDirection(Protein inProt){
			UnitVector lastVectorDirection = Param.ZERO_VECTOR;
			Coordinate last = inProt.getMap().get(inProt.getMap().size()-1);
			Coordinate nextToLast = inProt.getMap().get(inProt.getMap().size()-2);
			
			//get the component vectors
			//in the 2D square lattice one direction is 1 and the other 0
			int xVector = last.getX() - nextToLast.getX();
			int yVector = last.getY() - nextToLast.getY();
			
			//now assign the lastVectorDirection accordingly
			if(xVector == 0){//either north or south
				if(yVector < 0){//the previous residue was towards north
					lastVectorDirection =Param.j_VECTOR;
				}
				else if(yVector > 0){//the previous residue was towards south
					lastVectorDirection = Param.J_VECTOR;
				}
			}else if(yVector == 0){//either east or west
				if(xVector < 0){//the previous residue was towards east
					lastVectorDirection = Param.i_VECTOR;
				}
				else if(xVector > 0){//the previous residue was towards west
					lastVectorDirection = Param.I_VECTOR;
				}
			}else{ //there is an error because in the square lattice one must be 0
				System.out.println("Error in map coordinates for this protein");
				lastVectorDirection = Param.ZERO_VECTOR;
			}//end if-else to set the lastVectorDirection
			
			return lastVectorDirection;
		
		}//end getLastVectorDirection
	//----End getLastUnitVector-------------------------------------------------------
		
		
		
	//----findPossDirections-----------------------------------------------------------
		static public ArrayList<Direction> findPossDirections(Protein inProt, Param inParam){
			//the Coordinate of the last residue in the protein before adding 
			Coordinate lastCoord = inProt.getMap().get(inProt.getMap().size()-1);
			//make a temporary coordinate to use for testing
			Coordinate testCoord = new Coordinate();
			//set the last residue's direction
			UnitVector lastResDirection = getLastVectorDirection(inProt);
					
			//make a list of the possible directions
			ArrayList<Direction> possDirections = new ArrayList<Direction>();
			//array of all the impossible directions
			ArrayList<Direction> imposs = new ArrayList<Direction>();
			
			

			
			//determine which of the Directions the residue could go
			for(int i=0; i<inParam.getNumDirections(); i++){
				//reset the iterator for the map
				Iterator<Coordinate> mapIter = inProt.getMap().iterator();
				//test each possible direction by the last residues position relative to the one before it
				switch(Param.allDirections2D[i]){
					case L: switch(lastResDirection){
										case i:testCoord.setX(lastCoord.getX());
											   testCoord.setY(lastCoord.getY() - Param.LATTICE_UNIT);
											   //go throught the map coordinates
											   while(mapIter.hasNext() ){
												   //if any are identical to the proposed new one remove that from the possible directions
												  if(testCoord.equals(mapIter.next()) ){
													  //remove from possible directions
													  imposs.add(Param.LEFT);
													  break; //cut out of loop
												  }//end if
											   }//end while
											   
											break;
										case I:testCoord.setX(lastCoord.getX());
										       testCoord.setY(lastCoord.getY() + Param.LATTICE_UNIT);
											   while(mapIter.hasNext() ){
												  if(testCoord.equals(mapIter.next()) ){
													  imposs.add(Param.LEFT);
													  break;
												  }//end if
											   }//end while
											   
											break;
										case j:testCoord.setX(lastCoord.getX() + Param.LATTICE_UNIT);
										       testCoord.setY(lastCoord.getY());
											   while(mapIter.hasNext() ){
												  if(testCoord.equals(mapIter.next()) ){
													  imposs.add(Param.LEFT);
													  break;
												  }//end if
											   }//end while
											   
											break;
										case J:testCoord.setX(lastCoord.getX() - Param.LATTICE_UNIT);
										       testCoord.setY(lastCoord.getY());
											   while(mapIter.hasNext() ){
												  if(testCoord.equals(mapIter.next()) ){
													  imposs.add(Param.LEFT);
													  break;
												  }//end if
											   }//end while
											   
											break;
										default: System.out.println("Error in finding possible spots for elongation, null unit vector");
									    break;
										}//end internal Switch for L
						break; //case L
					case R: switch(lastResDirection){
										case i:testCoord.setX(lastCoord.getX());
										       testCoord.setY(lastCoord.getY() + Param.LATTICE_UNIT);
											   while(mapIter.hasNext() ){
												  if(testCoord.equals(mapIter.next()) ){
													  imposs.add(Param.RIGHT);
													  break;
												  }//end if
											   }//end while
											   
											break;
										case I:testCoord.setX(lastCoord.getX());
										       testCoord.setY(lastCoord.getY() - Param.LATTICE_UNIT);
											   while(mapIter.hasNext() ){
												  if(testCoord.equals(mapIter.next()) ){
													  imposs.add(Param.RIGHT);
													  break;
												  }//end if
											   }//end while
											   
											break;
										case j:testCoord.setX(lastCoord.getX() - Param.LATTICE_UNIT);
										       testCoord.setY(lastCoord.getY());
											   while(mapIter.hasNext() ){
												  if(testCoord.equals(mapIter.next()) ){
													  imposs.add(Param.RIGHT);
													  break;
												  }//end if
											   }//end while
											   
											break;
										case J:testCoord.setX(lastCoord.getX() + Param.LATTICE_UNIT);
										       testCoord.setY(lastCoord.getY());
											   while(mapIter.hasNext() ){
												  if(testCoord.equals(mapIter.next()) ){
													  imposs.add(Param.RIGHT);
													  break;
												  }//end if
											   }//end while
											   
											break;
										default: System.out.println("Error in finding possible spots for elongation, null unit vector");
									    break;
										}//end internal Switch for R
						break; //case R
					case S: switch(lastResDirection){
										case i:testCoord.setX(lastCoord.getX() - Param.LATTICE_UNIT);
										       testCoord.setY(lastCoord.getY());
											   while(mapIter.hasNext() ){
												  if(testCoord.equals(mapIter.next()) ){
													  imposs.add(Param.STRAIGHT);
													  break;
												  }//end if
											   }//end while
											   
											break;
										case I:testCoord.setX(lastCoord.getX() + Param.LATTICE_UNIT);
										       testCoord.setY(lastCoord.getY());
											   while(mapIter.hasNext() ){
												  if(testCoord.equals(mapIter.next()) ){
													  imposs.add(Param.STRAIGHT);
													  break;
												  }//end if
											   }//end while
											   
											break;
										case j:testCoord.setX(lastCoord.getX());
										       testCoord.setY(lastCoord.getY() - Param.LATTICE_UNIT);
											   while(mapIter.hasNext() ){
												  if(testCoord.equals(mapIter.next()) ){
													  imposs.add(Param.STRAIGHT);
													  break;
												  }//end if
											   }//end while
											   
											break;
										case J:testCoord.setX(lastCoord.getX());
										       testCoord.setY(lastCoord.getY() + Param.LATTICE_UNIT);
											   while(mapIter.hasNext() ){
												  if(testCoord.equals(mapIter.next()) ){
													  imposs.add(Param.STRAIGHT);
													  break;
												  }//end if
											   }//end while
											   
											break;
										default: System.out.println("Error in finding possible spots for elongation, null unit vector");
									    break;
										}//end internal Switch for S
						break; //case S
					default: System.out.println("Error in finding possible spots for elongation, null protein");
					    break;
											
					} //end main switch
				}//end for loop for each direction
			
			//add in all directions to possDirections
			for(int i=0; i < inParam.getNumDirections(); i++){
				possDirections.add(inParam.getAllDirections().get(i));
			}
					
			//now take the impossible ones out of possibleDirections
			//fill it with all Directions to remove the impossible ones later
			for(int i=0; i<imposs.size(); i++){
				for(int j=0; j < possDirections.size(); j++){
					if(imposs.get(i) == possDirections.get(j)){
						possDirections.remove(j);
					}//end if
				}//end inner for possDirections
			}//end outer for imposs

					
			return possDirections;
		}//end the find possible directions function
	//----End findPossDirections-------------------------------------------------------
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
	
//----makeBuildSeqs-----------------------------------------------------------

public static void makeBuildSeqs(ArrayList<Colony> expCols, Path expDir){
	
	//make build folders for selected builds		
			ListIterator<Colony> colIter = expCols.listIterator();
			Colony currCol = new Colony();
			int buildCount=1;
			while(colIter.hasNext()){
				currCol = colIter.next();
				ListIterator<ProteinBuild> buildIter = currCol.getBuilds().listIterator();
				while(buildIter.hasNext() && buildCount <= buildSeqNum){
					ProteinBuild currBuild = new ProteinBuild();
					currBuild = buildIter.next();
					try {
						Path buildDir;
						fillProtBuildImages(currBuild, currCol);
						//make a title for this build's directory
						String buildTitle = "Build_#" +buildCount +"_from_"+ currCol.getColTitle();
						buildDir = initBuildDir(expDir, buildTitle);
						
						ListIterator<BufferedImage> imgsIter = currBuild.getBuildStepsImages().listIterator();
						int imgCount = 1;
						while(imgsIter.hasNext()){
							
							ImageIO.write( imgsIter.next(),  "PNG", new File(buildDir.toString() + "\\step_" + imgCount + ".PNG") );
							imgCount++;
						}
						
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					buildCount++;
				}//end build iter
			}//end colony iter
	
}//end makeBuildSeqs
	
//----End makeBuildSeqs-------------------------------------------------------
	
//----fillProtBuildImages-----------------------------------------------------------
public static void fillProtBuildImages(ProteinBuild build, Colony colony) throws Exception{
	
	PNG_Footer footer = new PNG_Footer();
	
	for(int i=0; i < build.getBuildSteps().size(); i++){
		
		//set footer specifics
		footer.setSolutionStr("Solution rank:  " + i) ;
		footer.setEnergyStr("Energy:  " + build.getBuildSteps().get(i).getEnergy());
		footer.setSequenceInfoStr("Sequence: "+colony.getColParam().getSequenceName() + "  Length: "+ colony.getColParam().getProteinLength()+ "  Opt E:  "+ colony.getColParam().getOptimalEnergy());
		footer.setColInfoStr("Colony: " +colony.getColTitle());
		
		build.getBuildStepsImages().add( Protein.makeImage( build.getBuildSteps().get(i), footer ) );
		

}//end cycle through builds
	
}//End fillProtBuildImages
	
//----End fillProtBuildImages-------------------------------------------------------
	
//----initBuildDir-----------------------------------------------------------

public static Path initBuildDir(Path parentDir, String title){
	
	Path buildDir = Paths.get(parentDir + "\\" + title);
	System.out.println("Created: " + buildDir.toString());
	
	try{
	//create the directory
	Files.createDirectory(buildDir);
	
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
	
	
	return buildDir;
}//End initBuildDir
	
//----End initBuildDir-------------------------------------------------------

	
//----initExperiment-----------------------------------------------------------

	public static Path initExpDir(String expTitle) throws FileAlreadyExistsException{
		
		Path expDir = Paths.get(EXPERIMENT_BIN + expTitle );
		
		try {
			Files.createDirectory(expDir);
		} catch (FileAlreadyExistsException fe) {
			
			System.out.println("Experiment name \"" + expTitle + "\" already exists choose another");
			throw fe;	
			
		} catch (IOException e ) {
			e.printStackTrace();
		}
		
		return expDir;
	}//end initExperiment
//----End initExperiment-------------------------------------------------------
	
//----writeColonyDataFile-----------------------------------------------------------
	static public void writeColonyDataFile(Colony inColony, Path parentDir){ //colony to write results for and the parent directory path
		
		//create a colony directory
		Path colDir = Paths.get(parentDir + "\\" + inColony.getColTitle() );
		System.out.println("Created: " + colDir.toString());
		
		//the colony data file path
		Path colDataFile = Paths.get( colDir.toString() + "\\data.csv");
		System.out.println("Created: " + colDataFile.toString());
				
		try {
			//create the directory
			Files.createDirectory(colDir);
			//create the actual file
			Files.createFile(colDataFile);
			
			//write the file
			inColony.colResultFile(colDataFile);
			//create image files for all of the proteins
//DEBUG************************************************
			//inColony.makeProtImagePNGs(inColony.getSolutions(), colDir);
			
			
			//create image files for the solution if there is one
			if(inColony.isSolFound()){
				ImageIO.write(inColony.getSolImage(), "PNG", new File(colDir.toString() + "\\solutionImage.png") ); 
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
	}//end writeColonyDataFile
//----End writeColonyDataFile-------------------------------------------------------
	
//----writeTotalResults-----------------------------------------------------------
	public static void writeTotalResults(Path dir, ArrayList<Colony> cols){
		Path totalResults = Paths.get(dir.toString() +"\\total_results.csv");

		System.out.println("Created: " + totalResults.toString());
		
		Colony tempCol;
		
		
		Charset charset = Charset.forName("US-ASCII");
		try( BufferedWriter writer = Files.newBufferedWriter(totalResults, charset, WRITE, APPEND, CREATE) ){
			//the header line for each colony
			writer.write("wave,");
			
			//find the largest num waves
			int maxWaves=0;
			int tempWave = 0;
			ListIterator<Colony> colIter = cols.listIterator();
			while(colIter.hasNext()){
				tempWave = colIter.next().getWave();
				if(tempWave > maxWaves)
					maxWaves = tempWave;
			}//end search for max waves
			
			//write the waves on each column
			for(int i=0; i < maxWaves; i++){
				writer.write(" "+ (i+1) + ",");
			}
			writer.newLine();
			
			//write the wave averages
			writer.write("Wave Averages");
			writer.newLine();
			
			//now write each colonies name/row
			colIter = cols.listIterator();
			while(colIter.hasNext()){
				
				tempCol = colIter.next();
				writer.write(tempCol.getColTitle() + ",");
				
				//print out values
				ListIterator<Double> avgIter = tempCol.getWaveAverages().listIterator();
				while(avgIter.hasNext()){
					writer.write(" "+ avgIter.next().toString() +",");
				}
				writer.newLine();
				
			}//end write averages
			
			//write each colonies best solution energies
			writer.write("Solutions best to worst");
			writer.newLine();
			
			// write each colonies name/row
			colIter = cols.listIterator();
			while(colIter.hasNext()){
				
				tempCol = colIter.next();
				writer.write(tempCol.getColTitle() + ",");
				
				//print out values
				ListIterator<Protein> solIter = tempCol.getSolutions().listIterator();
				while(solIter.hasNext()){
					writer.write(" "+ solIter.next().getEnergy() +",");
				}
				writer.newLine();
				
			}//end write solutions energies in order
			
		}catch(IOException e){
			e.printStackTrace();
		}
	}// end writeTotalResults
//----End writeTotalResults-------------------------------------------------------
	
//----setExpTitleDateTime-----------------------------------------------------------
 	public static String setExpTitleDateTime(){
		//set the date and time and given exp name
		StringBuilder expDirNameBuild = new StringBuilder(100);
		expDirNameBuild.append("Exp_");
		String expDirName = "not init";
		//get the name of the experiment
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		try {
			System.out.println("Enter an experiment name: ");
			expDirNameBuild.append(in.readLine()) ;
			expDirNameBuild.append("_");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//date and time
		Calendar date = Calendar.getInstance();
		String timeStamp =  "Date_" + date.get(Calendar.DATE) +"."+ (date.get(Calendar.MONTH) +1) +"."+ date.get(Calendar.YEAR) +"_Time_"+ date.get(Calendar.HOUR_OF_DAY) +"."+ date.get(Calendar.MINUTE) +"."+
				date.get(Calendar.SECOND);
		//System.out.println(timeStamp);
		expDirNameBuild.append(timeStamp);
		
		expDirName = expDirNameBuild.toString();
		System.out.println(expDirName);
		return expDirName;
		
	}//end setExpTitleDateTime
//----End setExpTitleDateTime-------------------------------------------------------
	
}//end driver class
		


/*//OUTPUT CONSOLE
ListIterator<ArrayList<Protein>> waveIter = col1.getWaveBuilds().listIterator();

col1.printPheroMatrix("final pheromone matrix: ");

while(waveIter.hasNext()){
	ListIterator<Protein> antIter = waveIter.next().listIterator();
	while(antIter.hasNext()){
		Protein currProt = new Protein();
	
		currProt = antIter.next();
		System.out.println("Length: " + currProt.getLength());
		System.out.println("Energy: " + currProt.getEnergy());
		currProt.printCoordMap();
	}//end ant iterr
	
}//end wave iter
*/	