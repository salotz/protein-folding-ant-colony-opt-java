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

/* -------  Colony Class-------------------------------------------------------

	 
	
	ATTRIBUTES:
	
	
	CONSTRUCTORS:

	
	METHODS:
	GETTERS
	SETTERS
	
	
	
	
	Version 1.0 beta 
		
  ------------------------------------------------------------------------------
 */




package protein_folder;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import static java.nio.file.StandardOpenOption.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Random;

import javax.imageio.ImageIO;

import protein_folder.Param.*;

public class Colony {

//==================ATTRIBUTES======================================================
	//ATTRIBUTES
	private String colTitle;
	private int numAnts; //number of protein builds
	private ArrayList<ArrayList<Double>> pheroMatrix = new ArrayList<ArrayList<Double>>();  //the trails of pheromone
	private Param colParam;  //the set of paramters and sequence that this colony is working under/on
	private ArrayList<Protein> protBuilds = new ArrayList<Protein>(); //array of fully folded proteins
	private ArrayList<ArrayList<Protein>> waveBuilds = new ArrayList<ArrayList<Protein>>();
	private int wave;  //indicates the last wave or number to be performed
	
	private double selThresh;
	private BufferedImage solImage;
	private boolean solFound = false;
	
	//results attributes
	private ArrayList<Double> waveAverages = new ArrayList<Double>();
	private ArrayList<Protein> solutions = new ArrayList<Protein>();
	private ArrayList<ProteinBuild> builds = new ArrayList<ProteinBuild>();
	
//==================================================================================
	
	
	
//////__________Constructors_________________________________________________________

	
	//----no-arg constructor-----------------------------------------------------------

		public Colony() {
			// TODO Auto-generated constructor stub
	}
	//----end the no-arg constructor----------------------------------------------------
		
		
	//----Constructor-----------------------------------------------------------
	//arg constructor uses a Param object and produces waves of ants (inNumAnts = ants per wave) waves= number of waves
	public Colony(Param inParam, int inNumAnts, int waves, String inTitle){

		
		Protein solution =new Protein();
		
		//OPTIONS
		boolean random = true;
		boolean pheromoneTrails = false;
		boolean localSearch = false;
		boolean proteinBuild = false;
		
		colTitle = inTitle;
		colParam = inParam; //set the colony paramters to the input
		wave = 0;  //initialize the number of the current wave
		
		
		
		//turn on options
		for(Options inStratOpts: colParam.getOption()){
			switch(inStratOpts){
			case RANDOM: random = true;
				break;
			case PHEROMONE_TRAILS: pheromoneTrails=true;
				break;
			case LOCAL_SEARCH : localSearch = true;
				break;
			case PROTEIN_BUILD : proteinBuild = true;
				break;
			}
		}//end turning on options
		
		
		//if PHEROMONE_TRAILS option is on set to initPheromone
		if(pheromoneTrails){
			//set the selThresh attribute which will decrease with waves
			selThresh = colParam.getSelectionThreshold();
			///initialize the pheromone trails to tau= initPheromone from parameters
			for(int i=0; i< colParam.getProteinLength(); i++){
				ArrayList<Double> tauDirs = new ArrayList<Double>();
				for(int j=0; j< colParam.getNumDirections(); j++){
					tauDirs.add(colParam.getInitPheromone());
				}
				pheroMatrix.add(tauDirs);
			}//end initialize the pheromone trails
		
		}else {//set all tau to 1 which does not factor into equation
			
			///initialize the pheromone trails to tau= initPheromone from parameters
			for(int i=0; i< colParam.getProteinLength(); i++){
				ArrayList<Double> tauDirs = new ArrayList<Double>();
				for(int j=0; j< colParam.getNumDirections(); j++){
					tauDirs.add(Param.UNITY);
				}
				pheroMatrix.add(tauDirs);
			}//end initialize the pheromone trails
			
		}//end init pheromone trail values options guard
		
				
//DEBUG*********************************	
		printPheroMatrix("The initial pheroMatrix: ");
		
		
		
		
		
		//initialize waveBuilds
		for(int i=0; i < waves; i++){
			ArrayList<Protein> waveProtBuilds = new ArrayList<Protein>();
			waveBuilds.add(waveProtBuilds);
		}//end init waveBuilds
	
		
		//now make proteins in waves each with the input number of ants
		//wave iteration i
		int i=0;
		while(i < waves && solFound == false){
		
			System.out.println("Wave " + i + ":");
			//ant iterations j
			for(int j=0; j < inNumAnts; j++){
			
				
	//DEBUG********************************
				//System.out.println("Ant " + j + ":");
				
				//build a protein
				waveBuilds.get(i).add(antBuild(colParam.getResidueSeq()));
				
				//check to see if the recently made protein is a satisfactory solution
				if(waveBuilds.get(i).get(j).getEnergy() <= colParam.getSolutionCutoff()){
					
					solution = waveBuilds.get(i).get(j);
					solFound = true;
					
				}//solutionFound guard
			}//end ant loop for wave[i]
			
			//do the local seach if that option is on
			if(localSearch){
				//perform the LOCAL SEARCH on the proportion of ants that have been designated improving ants (i.e. do local search)
				for(int j=0; j < waveBuilds.get(i).size() * colParam.getProportionImprovingAnts(); j++){
					
					waveBuilds.get(i).set(j, localSearch( waveBuilds.get(i).get(j) ) );
					
					//check to see if the recently made protein is a satisfactory solution
					if(waveBuilds.get(i).get(j).getEnergy() <= colParam.getSolutionCutoff() ){
						solution = waveBuilds.get(i).get(j);
						solFound = true;
						
					}//solutionFound guard
					
				}//end local searches
			}//end local search option guard
				
				
			//only update if the option is on
			if(pheromoneTrails){
				//update the pheromone trails after each wave
				
				
	//DEBUG************************
				printPheroMatrix("pheromone before wave " + i + " : ");
				
				
				waveUpdatePheroMatrix();
				//reduce the selection threshold exponentially
				selThresh = Math.exp(-selThresh);
				if(selThresh < Param.MIN_SEL_THRESH)
					selThresh = Param.MIN_SEL_THRESH;
				
				
		//DEBUG************************
				printPheroMatrix("wave " + i +  " pheromone after: ");
				
			}//end PHEROMONE_TRAILS option guard
			
			i++;
		}//end waves loop
		
		if(solFound){
			System.out.println("An adequate solution was found");
			PNG_Footer solFooter = new PNG_Footer();
			//set footer specifics
			solFooter.setSolutionStr("Wave #  " + i) ;
			solFooter.setEnergyStr("Energy:  " + solution.getEnergy());
			solFooter.setSequenceInfoStr("Sequence: "+this.getColParam().getSequenceName() + "  Length: "+ this.getColParam().getProteinLength()+ "  Opt E:  "+ this.getColParam().getOptimalEnergy());
			solFooter.setColInfoStr("Colony: " +this.getColTitle());
			try {
				solImage = Protein.makeImage(solution, solFooter);
				
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		//calculate the statistics associated with the results
		this.calcStats();
		
		//create an ordered list of total proteins made in order of low to high energy
		//add them all to the ArrayList
		for(int m=0; m < waveBuilds.size(); m++){
			solutions.addAll(waveBuilds.get(m));
		}
		//sort them
		Collections.sort(solutions);
		
		
		
			
			

		
		
		
//DEBUG*********************************************
		//for(int i=0; i < solutions.size(); i++)
		//	System.out.println(solutions.get(i).getEnergy());
		
		
		
	}//end the Colony Constructor
//----End Constructor-------------------------------------------------------
	
//___________________________________________________________________________________
	
//==================METHODS======================================================
	
	

	public void makeProtImagePNGs(ArrayList<Protein> protList, Path filePath){
		ArrayList<BufferedImage> imgs = new ArrayList<BufferedImage>();
		int count=0;
		Protein currProt;
		try{
			
			ListIterator<Protein> solIter = protList.listIterator();
			count = 1;
			while(solIter.hasNext()){
				currProt = solIter.next();
				//create the footer for this protein
				PNG_Footer footer = new PNG_Footer();
				footer.setSolutionStr("Solution rank:  " + count) ;
				footer.setEnergyStr("Energy:  " + currProt.getEnergy());
				footer.setSequenceInfoStr("Sequence: "+colParam.getSequenceName() + "  Length: "+ colParam.getProteinLength()+ "  Opt E:  "+ colParam.getOptimalEnergy());
				footer.setColInfoStr("Colony: " +getColTitle());
				
				imgs.add( Protein.makeImage(currProt, footer) );
				
				count++;
			}
			
			ListIterator<BufferedImage> imgsIter = imgs.listIterator();
			count = 1;
			while(imgsIter.hasNext()){
				
				ImageIO.write( imgsIter.next(),  "PNG", new File(filePath.toString() + "\\solution_" + count + ".PNG") );
				count++;
			}
			
			
			
		}catch (IOException ie) {
			ie.printStackTrace();
		}
		catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}//end the making images try-catch
		
		
	}//end makeProtImagePNGs
	
//----End makeProtImagePNGs-------------------------------------------------------
	

//----localSearch-----------------------------------------------------------

	//takes as an argument the proteins current conformation
	//returns either the same protein conformation or one with an improved energy conformation
	//will attempt to make better conformations either until it finds a better one
	//of the localSearchMaxIter number has been reached and none have been found
	//the probability that a residue changes its relative next direction is given in the colParam
	public Protein localSearch(Protein currConf){
		
		Protein outProt = new Protein(); //the protein the function will return
		Random rand = new Random();
		int mutResNum = 0;  // the residue index that the local search long range movement starts at
		boolean feasible = false; //flag raised when the current residues direction is feasible
		Direction tempDir = Param.END; 
		boolean abort = false; //raised when the current local search is problematic and should be aborted
		boolean solution = false; //raised when a conformation is found that is better than the current one
		
		int i=0;
		//create new conformations until a better one is found or the limit for searches is up
		while(i < colParam.getLocalSearchMaxIter()  &&   !solution){
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
				feasibleNewDirs.addAll( findPossDirections(newConf) );
				
				//check to see if the currConf previous position is infeasible
				feasible = false; //if a match is found in possDirs then it is
				for(int k=0; k < feasibleNewDirs.size(); k++){
					if(currConf.getSequence().get(mutResNum + j).getNextResDirection() == feasibleNewDirs.get(k)){
						feasible = true;
					}
				}//end check for feasibility
				
				//if probability wills it or the previous direction is infeasible make new direction
				if( rand.nextDouble() <= colParam.getProbKeepConfLocalSearch() ||  !feasible) {
					
					//set the proposed fold to tempDirection
					tempDir = addFoldLocalSearch(currConf.getSequence().get(mutResNum + j).getType(), newConf);
					
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
	
	
//----initColDir-----------------------------------------------------------

	public static Path initColDir(Path expDir, String colTitle) throws FileAlreadyExistsException{
		
		Path colDir = Paths.get(expDir + colTitle );
		
		try {
			Files.createDirectory(colDir);
		} catch (FileAlreadyExistsException fe) {
			
			System.out.println("Experiment name \"" + colTitle + "\" already exists choose another");
			throw fe;	
			
		} catch (IOException e ) {
			e.printStackTrace();
		}
		
		return colDir;
	}//end initColDir
//----End initColDir-------------------------------------------------------
	
//----colResultFile-----------------------------------------------------------	
	public void colResultFile(Path outputFile){
		Charset charset = Charset.forName("US-ASCII");
	
		try( BufferedWriter writer = Files.newBufferedWriter(outputFile, charset, WRITE, CREATE, APPEND) ){
			
			//print the colony parameter identifier
			writer.write(colParam.getTitle());
			writer.newLine();
			
			//Print the cell header line
			writer.write("Wave, Build, Energy, Average, Std Dev, Sequence Name, Length, Optimal Energy, Strategy Options, Number of Coordinate Directions," +
					"Inverse Temperature, alpha, beta, rho, theta, init tau, Selection Threshold, Probability a Feasible Conformation is kept in Local Search, Proportion Improving Ants");
			writer.newLine();
			
			
			//print the values of parameters for this colony	
			writer.write("  ,  ,  ,  ,  ,"+ colParam.getSequenceName() +","+ colParam.getProteinLength() +","+ colParam.getOptimalEnergy() +","+ colParam.getOption() +","+
						colParam.getNumDirections() +","+ colParam.getInvTemp() +","+ colParam.getAlpha() +","+ colParam.getBeta() +","+ colParam.getRho() +","+ colParam.getTheta() +","+ 
					colParam.getInitPheromone() +","+ colParam.getSelectionThreshold() +","+ colParam.getProbKeepConfLocalSearch() +","+ colParam.getProportionImprovingAnts()
				);
			writer.newLine();
			
			ListIterator<ArrayList<Protein>> waveIter = waveBuilds.listIterator();			
			int waveCount=1;
			int antCount=1;
						
			while(waveIter.hasNext()){
				ListIterator<Protein> antIter = waveIter.next().listIterator();
				writer.write(waveCount +",");
				writer.newLine();
				antCount=0;
				while(antIter.hasNext()){
					Protein currProt = new Protein();
					currProt = antIter.next();
					writer.write(" ," + antCount +","+ currProt.getEnergy());
					writer.newLine();
					antCount++;
				}//end ant iterr
				
				//write the average
				writer.write(" , , ," + waveAverages.get(waveCount-1) +",");
				writer.newLine();
								
				waveCount++;
			}//end wave iter
			
	} catch (IOException e){
			e.printStackTrace();
		}
		
	}//end colResultFile
//----End colResultFile-------------------------------------------------------
	
//----calcStats-----------------------------------------------------------
	public void calcStats(){
		
		ListIterator<ArrayList<Protein>> waveIter = waveBuilds.listIterator();			
		int antCount=1;
		double waveEnergy = 0;
		double waveAvg = 0;
		
		//find the average energies for each wave
		//find the best solution from all waves
		while(waveIter.hasNext()){
			ListIterator<Protein> antIter = waveIter.next().listIterator();
			antCount=0;
			while(antIter.hasNext()){
				Protein currProt = new Protein();
				currProt = antIter.next();
						
				//the sum of energy for the wave
				waveEnergy += currProt.getEnergy();
				antCount++;
			}//end ant iterr
			
			//calculate average energy for the wave
			waveAvg = waveEnergy/antCount;
			//store it in the results array for the colony
			waveAverages.add(waveAvg);
			//calculate standard deviation for the wave
			//STUB
			
			//updates
			waveAvg = 0;
			waveEnergy = 0;
			
		}//end wave Iter
	}//end calcStats
//----End calcStats-------------------------------------------------------
	
	
	
	
	
//----waveUpdatePheroMatrix-----------------------------------------------------------	
	//this function updates the pheromone trails after each wave
	public void waveUpdatePheroMatrix(){
		Random rand = new Random();
		ArrayList<Integer> energies = new ArrayList<Integer>();//to store the energies for the last wave of builds
		int cutoff=0;  //used to find the energy cutoff for allowing a protein build to update the pheromone matrix
		int cutoffE=0;  //the energy a protein must have to be allowed to update the pheroMatrix
		Protein tempProt;
		double tempTau = 0;
		ArrayList<Double> minTaus = new ArrayList<Double>();
		ArrayList<Double> sumTaus = new ArrayList<Double>();
		ArrayList<Double> normMinTaus = new ArrayList<Double>();
		
		ArrayList<Integer> normList = new ArrayList<Integer>();
		int maxIndex = 0;
		double maxTau = 0;
		ArrayList<Integer> maxTieList = new ArrayList<Integer>();
		
		//initialize minTaus and sumTaus
		for(int i=0; i< pheroMatrix.size(); i++){
			minTaus.add(Param.THOUSAND); //set all min taus to 1000 to ensure that all actual tau values are less than this
			sumTaus.add(Param.ZERO);
		}//end for
		
		
//DEBUG**************
		printPheroMatrix("before evaporation: ");
		
		
		//First: "evaporate" pheromone by applying the persistence (rho) multiplier to each tau
		for(int i=0; i< pheroMatrix.size(); i++){
			for(int j=0; j< pheroMatrix.get(i).size(); j++){
				//set it to 
				pheroMatrix.get(i).set(j, colParam.getRho() * pheroMatrix.get(i).get(j) );
				
			}//end inner for
		}//end outer for
		
		
//DEBUG**************
		printPheroMatrix("after evaporation: ");
		
		
		//Second: we must find the low energy conformations using the selection threshold value 
		//using the last wave
		//load the energies array list it with the energies of the last waves builds
		for(int i=0; i < waveBuilds.get(wave).size(); i++){
			energies.add(waveBuilds.get(wave).get(i).getEnergy());
		}//end load energies
		
		//now sort it lowest to highest energy
		 Collections.sort(energies);
		 
		 //set the cutoff value for the builds with 10% lowest energies
		 cutoff = (int) (energies.size()*colParam.getSelectionThreshold());
		 //automatically set cutoff to 1 if it is 0
		 if(cutoff == 0)
			 cutoff = 1;
		 cutoffE = energies.get(cutoff);
				 
		 
//DEBUG*********************************
			printPheroMatrix("Before significant protein: ");
		 
		 
		//now update the pheroMatrix based on these selected protein builds
			for(int i=0; i < waveBuilds.get(wave).size(); i++){
				//if a protein build has energy < cutoff energy then use it to update the matrix
				tempProt = waveBuilds.get(wave).get(i);
				if(tempProt.getEnergy() <= cutoffE){
					//use the single protein updater
					updatePheroMatrix(tempProt);
					
//DEBUG*********************************
					printPheroMatrix("After significant protein " +i+ " : ");
						
				}//end conditions for updateing the pheroMatrix
			}//end cycle through last builds for
			
//DEBUG*********************************
				printPheroMatrix("After significant proteins: ");
				
				
			//Third: renormalize the pheroMatrix so that sumTau for each seqNum = 1
			

				//A: find sumTau for each seqNum
				for(int i=0; i < pheroMatrix.size(); i++){					
					for(int j =0; j < pheroMatrix.get(i).size(); j++){
						
						tempTau = pheroMatrix.get(i).get(j);
						//add this tau to the sumTau
						sumTaus.set(i, sumTaus.get(i) + tempTau);
						
					}//end inner j loop	
					
				}//end find min loop outer i
				
				//B: normalize every seqNum tau to each seqNum sumTau
	//DEBUG*********************************
				printPheroMatrix("Before normalization :" );
				
				for(int i=0; i < pheroMatrix.size(); i++){
					for(int j =0; j < pheroMatrix.get(i).size(); j++){
						
						//reset each tau to a normalized value
						pheroMatrix.get(i).set(j, pheroMatrix.get(i).get(j) / sumTaus.get(i) ); 
						
					}//end inner j loop		
										
				}//end find min loop outer i
		//DEBUG*********************************
				printPheroMatrix("after normalization :" );
				
			//FOURTH: Implement the MIN-MAX adaptation system if any tau < THETA adapt extreme values
				
				//A: find the minTaus for each seqNum
				for(int i=0; i < pheroMatrix.size(); i++){
										
					for(int j =0; j < pheroMatrix.get(i).size(); j++){
						
						tempTau = pheroMatrix.get(i).get(j);
						
						//if the current tau is less than the minTau set minTau to that
						if(tempTau < minTaus.get(i)){
							minTaus.set(i, tempTau);
						}//end if for minTau
						
					}//end inner j loop			
				}//end find min loop outer i
				//end 4A
				
				
				//B : Go through and for each seq number with minTau < THETA, adapt extreme values
				//DEBUG*********************************
				printPheroMatrix("Before adaptation :" );
				for(int i=0; i < pheroMatrix.size(); i++){
					
					//if the minTau has dropped below THETA then we need to normalize 
					if( minTaus.get(i) < colParam.getTheta() ){
						
						

						
						
						//re-initialize maxTau
						maxTau = 0;
						//find the maxTau value, check for ties
						for(int j=0; j < pheroMatrix.get(i).size(); j++){
							//set the tempTau
							tempTau = pheroMatrix.get(i).get(j);
							//if the current tempTau is greater than the maxTau set MaxTau to that
							if(tempTau > maxTau){
								
								maxTau = tempTau;
								//for ties
								maxIndex = j; //set the index of the maxTau so we can record it later in the case of a tie
								//the old tie is not valid as highest anymore so clear it
								maxTieList.clear();
								
							}//end if for maxTau
							//if there is a tie for max tau record which ones were tied 
							//(i.e. 1,2, and 3 corresponding later to directions say S, L, and R)
							else if(tempTau == maxTau  && j != 0){ //maxIndex != 0 to protect ourselves
								//include the first tied max value in the list unless it was already incorporated
								if(!maxTieList.contains(maxIndex))
									maxTieList.add(maxIndex);
								//then add the second tied value in the list
								maxTieList.add(j);
								
							}//end tie for maxTau
							
						}//end find maxTau or tie list
						
						//decide which tied direction to choose to change uniformly randomly
						if(!maxTieList.isEmpty()){
							maxIndex = Math.abs( rand.nextInt() % maxTieList.size() );
						}//end if for choosing tied values
						
						//decrease maxTau := maxTau - (THETA-minTau)
						try{
							
						pheroMatrix.get(i).set(maxIndex, maxTau - Math.abs( colParam.getTheta() - minTaus.get(i) ) );
						}
						catch(ArrayIndexOutOfBoundsException e){
							e.printStackTrace();
							
						}
						
						// then set minTau to THETA  for tied minTaus all are set to THETA
						for(int j=0; j < pheroMatrix.get(i).size(); j++){
							tempTau = pheroMatrix.get(i).get(j);
							
							//so if any tau < THETA && tau == minTau set tau=THETA
							if(tempTau < colParam.getTheta() && tempTau == minTaus.get(i)){
								
								pheroMatrix.get(i).set(j, colParam.getTheta());
								
							}//end check for tied minTau values
						}//end set all min taus to theta
						
	
						
					}//end check for normMinTaus < theta
					
			
					
				}//end renormalize loop
		
				//DEBUG*********************************
				printPheroMatrix("after adaptation :" );
				
	}//end update pheroMatrix for wave
//----End waveUpdatePheroMatrix-------------------------------------------------------
	

	
//----updatePheroMatrix-----------------------------------------------------------
	//updatePheroMatrix uses a given protein to increase pheromone
	public void updatePheroMatrix(Protein inProt){
		
		//calculate the relative solution quality
		Double relSolutionQuality = inProt.getEnergy() / colParam.getOptimalEnergy();
		//iterator for the sequences
		ListIterator<Residue> seqIter = inProt.getSequence().listIterator();
		//initialize the next direction by getting the first direction
		int i = 0;
		
		//initialize entrance to while loop
		Direction dir = seqIter.next().getNextResDirection();
		//check to makes sure protein is there
		if(dir == Param.END){
			System.out.println("Error empty protein in updatePheroMatrix");
		}//end check for non-empty protein
		// find the index of that direction on the allDirections2D constant array
		//use that to update the corresponding index position in pheroMatrix
		for(int j=0; j <  colParam.getAllDirections().size(); j++){
			if(dir == colParam.getAllDirections().get(j)){
				
				//set newTau = tau + relSolutionQuality
				pheroMatrix.get(i).set(j, pheroMatrix.get(i).get(j) + relSolutionQuality);
				
			}//end get index of direction
		}//end loop through directions
		i++;
		
		
		
		
		//Go through each seq number in the pheromone trail
		while(seqIter.hasNext() && dir != Param.END){
			
			//get the direction of the next residue
			dir = seqIter.next().getNextResDirection();
			
			//now find the index of that direction on the allDirections2D constant array
			//use that to update the corresponding index position in pheroMatrix
			for(int j=0; j <  colParam.getAllDirections().size(); j++){
				if(dir == colParam.getAllDirections().get(j)){
					
					//set newTau = tau + relSolutionQuality
					pheroMatrix.get(i).set(j, pheroMatrix.get(i).get(j) + relSolutionQuality);
					
				}//end get index of direction
			}//end loop through directions	
			i++;
		}//end out for i and seqIter
		
	}//end update pheroMatrix for a single protein build
//----End updatePheroMatrix-------------------------------------------------------
	
	
	
//----printPheroMatrix-----------------------------------------------------------
	public void printPheroMatrix(String title){
		System.out.println(title); //print the title
		
		for(int i=0; i < pheroMatrix.size(); i++){
			System.out.print("{ ");
			for(int j=0; j< pheroMatrix.get(i).size(); j++){
				System.out.print(pheroMatrix.get(i).get(j) + ", ");
			}
			System.out.println("}");
		}//end seq num loop for
	}//end print pheroMatrix
//----End printPheroMatrix-------------------------------------------------------

	
	
//----antBuild-----------------------------------------------------------
	//construct an entire protein
	//sequence size must be >2 to avoid errors
	public Protein antBuild(ArrayList<ResidueType> inSeq){
		//iterator for residue sequence list
		ListIterator<ResidueType> seqIter = inSeq.listIterator();
		
		//check to make sure the sequence is non-empty
		if(inSeq.size() > 2){
			Protein newProt = new Protein(seqIter.next()); //start the protein with it's first residue
			newProt.addResidue(Param.STRAIGHT, seqIter.next()); //then the second the direction does not matter
			
			//add build to list
			/*ProteinBuild newBuild = new ProteinBuild();
			builds.add(newBuild);
			//make first buildStep the two residue starting sequence
			newBuild.getBuildSteps().add(new Protein(newProt) );*/
			
			//now go through the rest
			while(seqIter.hasNext()){
				ResidueType tempType = seqIter.next();				
				//a variable for the next residues direction
				Direction newDirection = Param.END;
				
				//set the new direction by calling the function that can find possible folds and uses probabilistic calculations
				newDirection = addFold(tempType, newProt);
				//if it is the null direction E than we must do attrition
				if(newDirection == Param.END){
					
					//do attrition and reset the newProt
					newProt = doAttrition(newProt);
					
					//add new attrition protein to the build list
					//newBuild.getBuildSteps().add(new Protein(newProt) );
					
					//check to make sure attrition worked
					//will return protein of length 0 if it did not
					if(newProt.getLength() == 0){
						
		//DEBUG***********************************************
					//	System.out.println("FATAL ATTRITION");
						
						//restart the protein and iterator
						seqIter = inSeq.listIterator();
						newProt = new Protein(seqIter.next()); //first residue
						newProt.addResidue(Param.STRAIGHT, seqIter.next());  //second residue
						//add new attrition protein to the build list
						//newBuild.getBuildSteps().add(new Protein(newProt) );
						
					}//end check for successful attrition
					
					//if attrition was done move the list Iterator back one so it can redo that residue
					seqIter.previous();
					
				}else{ //it is fine and a new residue can be added
					newProt.addResidue(newDirection, tempType); //add a residue to this protein build
					
					//add this step to the build object list
					//newBuild.getBuildSteps().add(new Protein(newProt) );
					
				}//end if else for doing attrition or adding a residue
				
				
			}//end iteration while
				return newProt;
		}else{
			Protein newProt = new Protein(Param.NULL);
			System.out.println("The input sequence to this ant is less than 2, please choose a longer sequence");
			
			return newProt;
		}//end check for existence of residueSequence
		
		
		
	}//end ant build function	
//----End antBuild-------------------------------------------------------
	
	
	
//----addFold-----------------------------------------------------------	
	//adds a fold to the protein by finding the options, making a choice, and then returning the direction
	// this is a large function with multiple phases which I have included all in one for now because the sub-functions
	//require a lot of parameters and it might be easier to house them under one function
	//Phase 1: find all the spots that the new residue could possibly be put
	//Phase 2: if no spots do attrition
	//Phase 3: Make a probabilistic choice based on the pheromone and thermo favorability
	//returns the direction of the next residue
	public Direction addFold(ResidueType inType, Protein inProt){
		Random rand = new Random();
		//the direction we will return
		Direction newDirection = Param.END;
		double choiceNum = rand.nextDouble(); //used to pick a direction (0.0 <= choiceNum < 1.0)
		double cumProb = 0; //used for calculating cumulative probabilities
		
		//make a list of the possible directions and fill it
		ArrayList<Direction> possDirections = new ArrayList<Direction>();
		possDirections.addAll(findPossDirections(inProt));
		
		//the cumulative distribution boundaries, first element is 0
		Double[] probs = new Double[possDirections.size() + 1];
		probs[0] = cumProb; //because it is 0 right now
		
		
		//if there are no possible directions we must do attrition in the built protein
		//so return the null direction to signal to the ant to do attrition
		if( possDirections.isEmpty()){
			
//DEBUG***********************************************			
			//System.out.println("no possible directions attrition triggered");
			
			
			
			return Param.END;
		}
		else{// we choose which direction to return based on the thermodynamic model
		
			//first create the cumulative probability distribution
			cumProb=0;
			for(int i=0; i < possDirections.size(); i++){
				probs[i+1] = (cumProb + calcFoldProb(inProt, inType, possDirections.get(i)));  //DEBUG take out the first parameter
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
	}//end addfold method
//----End addFold-------------------------------------------------------

//----addFoldLocalSearch-----------------------------------------------------------	
	//adds a fold to the protein by finding the options, making a choice, and then returning the direction
	// this is a large function with multiple phases which I have included all in one for now because the sub-functions
	//require a lot of parameters and it might be easier to house them under one function
	//Phase 1: find all the spots that the new residue could possibly be put
	//Phase 2: if no spots do attrition
	//Phase 3: Make a probabilistic choice based on the pheromone and thermo favorability
	//returns the direction of the next residue
	public Direction addFoldLocalSearch(ResidueType inType, Protein inProt){
		Random rand = new Random();
		//the direction we will return
		Direction newDirection = Param.END;
		double choiceNum = rand.nextDouble(); //used to pick a direction (0.0 <= choiceNum < 1.0)
		double cumProb = 0; //used for calculating cumulative probabilities
		
		//make a list of the possible directions and fill it
		ArrayList<Direction> possDirections = new ArrayList<Direction>();
		possDirections.addAll(findPossDirections(inProt));
		
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
				probs[i+1] = (cumProb + calcFoldProbLocalSearch(inProt, inType, possDirections.get(i)));  //DEBUG take out the first parameter
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
	
//----calcFoldProb-----------------------------------------------------------
	//the calculate fold probability method
	public double calcFoldProb(Protein inProt, ResidueType inType, Direction foldDirection){
		double outProb = 0.0;
		double tau = 0.0; //pheromone level
		double eta = 0;
		int newHH = 0; //the new H-H bonds created by the input fold
		double weightNumerator = 0.0; //the weight of the fold of interest i.e. the numerator of the probability function
		double weightTotal = 0.0;  //the denominator in the probability function
		int tauDirIndex = 0;
		
		
		
		//first find the numerator or the value of the fold in question
		
		
		//find the tau value at the proposed fold
			//get the pheroMatrix inner index (seq position next direction array)
			for(int i=0; i < Param.allDirections2D.length; i++){
				if(foldDirection == Param.allDirections2D[i]){
					tauDirIndex = i;
				}
			}//end tauDir index get
		tau = pheroMatrix.get(inProt.getLength() -1 ).get(tauDirIndex);
		
		
		//protein that is constructed to calculate the energy from the input fold
		Protein testProtNum = new Protein(inProt);
		//array for possible folds
		ArrayList<Direction> possFolds = new ArrayList<Direction>();
		
		//calculate the number of new H-H bonds by making the test protein and calculating the energy difference
		testProtNum.addResidue(foldDirection, inType);
		// the getEnergy automatically updates the energy -- testProt.updateEnergy();
		newHH = testProtNum.getEnergy() - inProt.getEnergy();
		//calculate the eta value
		eta = Math.exp(-colParam.getInvTemp() * newHH);
		//find the tau value from the pheroMatrix

		//save to the numerator value
		weightNumerator = Math.pow(tau, colParam.getAlpha())*Math.pow(eta, colParam.getBeta());
		
	//DEBUG*********************************************
		//System.out.println("The numerator : " + weightNumerator);
		
		
		//NEXT SUM OF ALL POSSIBLE FOLDS
		//calculate the sum of all the possible folds
		possFolds = findPossDirections(inProt);
		//calculate the weights and add them up
		for(int i=0; i<possFolds.size(); i++){
			
			//get the tau value for the move
			//get the pheroMatrix inner index (seq position next direction array)
			for(int j=0; j < Param.allDirections2D.length; j++){
				if(possFolds.get(i) == Param.allDirections2D[j]){
					tauDirIndex = j;
				}
			}//end tauDir index get
			tau = pheroMatrix.get(inProt.getLength() -1 ).get(tauDirIndex);
			
			
			
			Protein testProtTemp = new Protein(inProt);
			//calculate the number of new H-H bonds by making the test protein and calculating the energy difference
			testProtTemp.addResidue(possFolds.get(i), inType);
			// the getEnergy automatically updates the energy -- testProt.updateEnergy();
			newHH = testProtTemp.getEnergy() - inProt.getEnergy();
			//calculate the eta value
			eta = Math.exp(-colParam.getInvTemp() * newHH);
			//find the tau value
			//STUBBED OUT
			
			//add it to the weight total
			weightTotal = weightTotal + Math.pow(tau, colParam.getAlpha())*Math.pow(eta, colParam.getBeta());
			
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
	}//end calcFoldProb
//----End calcFoldProb-------------------------------------------------------
	
//----calcFoldProbLocalSearch-----------------------------------------------------------
	//the calculate fold probability method
	public double calcFoldProbLocalSearch(Protein inProt, ResidueType inType, Direction foldDirection){
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
		eta = Math.exp(-colParam.getInvTemp() * newHH);
		//find the tau value from the pheroMatrix

		//save to the numerator value
		weightNumerator = eta;
		
	//DEBUG*********************************************
		//System.out.println("The numerator : " + weightNumerator);
		
		
		//NEXT SUM OF ALL POSSIBLE FOLDS
		//calculate the sum of all the possible folds
		possFolds = findPossDirections(inProt);
		//calculate the weights and add them up
		for(int i=0; i<possFolds.size(); i++){
							
			Protein testProtTemp = new Protein(inProt);
			//calculate the number of new H-H bonds by making the test protein and calculating the energy difference
			testProtTemp.addResidue(possFolds.get(i), inType);
			// the getEnergy automatically updates the energy -- testProt.updateEnergy();
			newHH = testProtTemp.getEnergy() - inProt.getEnergy();
			//calculate the eta value
			eta = Math.exp(-colParam.getInvTemp() * newHH);
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
	
//----doAttrition-----------------------------------------------------------
	//attrition function
	//this function is included in the ant colony class because I want to keep the protein class
	//basic in form so that it could be used in other building methods
	//attrition is part of the ant colony specific folding program
	//takes in a protein as an argument and returns the protein with attrition done to it
	public Protein doAttrition(Protein inProt){
		Random rand = new Random();
		//create a copy of inProt to modify
		Protein outProt = new Protein(inProt); 
		
		int lostRes = 0;
		int resIndex = 0;
		ArrayList<Direction> possDirections = new ArrayList<Direction>();
		
		//remove half of the protein
		for(int i=0; i < (inProt.getLength()) / 2; i++){
			outProt.removeResidue();
		}//end for removing residues
		
		//now replace the lost residues with randomly placed residues
		lostRes = inProt.getLength()- outProt.getLength(); //the number of residues to be replaced
		resIndex = outProt.getLength();  //the index of the residue in the sequence to be replaced next
		
		for(int i=0; i< lostRes; i++){
			//first look for possible places to go
			possDirections.clear(); //clear it from the old possible moves
			possDirections.addAll(findPossDirections(outProt));	
			
			//there may be a need for attrition again so recursively call it if there are no possible directions again
			if(possDirections.isEmpty()){
				outProt = doAttrition(outProt);
				possDirections.clear(); //clear it from the old possible moves
				possDirections.addAll(findPossDirections(outProt));	
				
				
//DEBUG***********************************************
				//System.out.println("successful recursive attrition");
				
				
			}
		
			//if there are still no directions to go the protein needs to be rebuilt from the beginning
			if(possDirections.isEmpty()){
				outProt.setLength(0);
				break;
			}
			//if not it is ok to continue
			else{
			//then pick one direction at random and place the next residueType in the residueSeq
			outProt.addResidue(possDirections.get(Math.abs(rand.nextInt() % possDirections.size())), colParam.getResidueSeq().get(resIndex));
			}
			//increase the resIndex to the next residue in residueSeq
			resIndex++;
		}//end replace residues
		
		return outProt;
	}//end attrition function
//----End doAttrition-------------------------------------------------------
	
	
	
//----getLastUnitVector-----------------------------------------------------------
	//get the unit vector direction to the last residue
	public UnitVector getLastVectorDirection(Protein inProt){
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
	public ArrayList<Direction> findPossDirections(Protein inProt){
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
		for(int i=0; i<colParam.getNumDirections(); i++){
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
		for(int i=0; i < colParam.getNumDirections(); i++){
			possDirections.add(colParam.getAllDirections().get(i));
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
	
	
	
//////++++++++GETTERS & SETTERS+++++++++++++++++++++++++++++++++++++++++++++++++++++

	public boolean isSolFound() {
		return solFound;
	}

	public void setSolFound(boolean solFound) {
		this.solFound = solFound;
	}

	public double getSelThresh() {
		return selThresh;
	}

	public void setSelThresh(double selThresh) {
		this.selThresh = selThresh;
	}

	public BufferedImage getSolImage() {
		return solImage;
	}

	public void setSolImage(BufferedImage solImage) {
		this.solImage = solImage;
	}
	public ArrayList<ProteinBuild> getBuilds() {
		return builds;
	}

	public void setBuilds(ArrayList<ProteinBuild> builds) {
		this.builds = builds;
	}
	public ArrayList<Protein> getSolutions() {
		return solutions;
	}
	public void setSolutions(ArrayList<Protein> solutions) {
		this.solutions = solutions;
	}
	public ArrayList<Double> getWaveAverages() {
		return waveAverages;
	}
	public void setWaveAverages(ArrayList<Double> waveAverages) {
		this.waveAverages = waveAverages;
	}
	public String getColTitle() {
		return colTitle;
	}
	public void setColTitle(String colTitle) {
		this.colTitle = colTitle;
	}
	public int getNumAnts() {
		return numAnts;
	}
	public void setNumAnts(int numAnts) {
		this.numAnts = numAnts;
	}
	public ArrayList<ArrayList<Double>> getPheroMatrix() {
		return pheroMatrix;
	}
	public void setPheroMatrix(ArrayList<ArrayList<Double>> pheroMatrix) {
		this.pheroMatrix = pheroMatrix;
	}
	public ArrayList<Protein> getProtBuilds() {
		return protBuilds;
	}
	public void setProtBuilds(ArrayList<Protein> protBuilds) {
		this.protBuilds = protBuilds;
	}
	public Param getColParam() {
		return colParam;
	}
	public void setColParam(Param colParam) {
		this.colParam = colParam;
	}
	public ArrayList<ArrayList<Protein>> getWaveBuilds() {
		return waveBuilds;
	}
	public void setWaveBuilds(ArrayList<ArrayList<Protein>> waveBuilds) {
		this.waveBuilds = waveBuilds;
	}
	public int getWave() {
		return wave;
	}
	public void setWave(int wave) {
		this.wave = wave;
	}

//////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
}//end colony class

//OLD Stubs
//METHODS
	//No-arg Constructor
	/*public Colony(){
		numAnts = 1;
		colParam = new Param();

		
		
		//create as many protein builds as dictated by numAnts
		//and put them into the protBuilds arrayList
		for(int i=0; i<numAnts; i++){
			
			//fill the protBuilds with builds from the next ant
			protBuilds.add(antBuild(colParam.getResidueSeq()));
			
	//DEBUG
			System.out.println("Protein made");
			
			
			//update the pheromone matrix
			//STUBOUT
		}//end for
	
	}//end no-arg constructor
	*/
	
	
	//default parameters argumented constructor for number of ants and sequence (from Param)
	/*	
	public Colony(int seqNum, int inNumAnts){
		numAnts = inNumAnts;
		//create the param object
		colParam = new Param(seqNum);
		
		///initialize the pheromone trails to tau=1.0
		for(int i=0; i< colParam.getProteinLength(); i++){
			ArrayList<Double> tauDirs = new ArrayList<Double>();
			for(int j=0; j< colParam.getNumDirections(); j++){
				tauDirs.add(Param.UNITY);
			}
			pheroMatrix.add(tauDirs);
		}//end initialize the pheromone trails
		
		//Now actually build the proteins
		for(int i=0; i < numAnts; i++){
			//make a protein
			protBuilds.add(antBuild(colParam.getResidueSeq()));
			
			//update the pheromones
			
	//STUBBED OUT    pheroMatrix.update(protBuilds.get(i));
			
		}//end make antBuilds loop
		
	}//end default argumented constructor
	*/
	
	//an argumented constructor that uses an already constructed Param object for it's argument
	/*
	public Colony(Param inParam, int inNumAnts){
		numAnts = inNumAnts; //set the number of ants to use total
		colParam = inParam; //set the param object
		
		///initialize the pheromone trails to tau=initPheromone
				for(int i=0; i< colParam.getProteinLength(); i++){
					ArrayList<Double> tauDirs = new ArrayList<Double>();
					for(int j=0; j< colParam.getNumDirections(); j++){
						tauDirs.add(colParam.getInitPheromone());
					}
					pheroMatrix.add(tauDirs);
				}//end initialize the pheromone trails
		
		
		for(int i=0; i < numAnts; i++){
			//make a protein
			protBuilds.add(antBuild(colParam.getResidueSeq()));
			
			//update the pheromones
			
	//STUBBED OUT    pheroMatrix.update(protBuilds.get(i));
			//***********************************
			
		}//end make antBuilds loop
		
	}//end Param specified constructor
*/
