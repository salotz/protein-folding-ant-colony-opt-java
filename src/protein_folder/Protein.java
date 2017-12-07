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

/* ------- Protein Class-------------------------------------------------------

	The protein class is intended to be a modular class that different folding 
	strategies can be applied to. 
	
	ATTRIBUTES:
	length -- merely holds the number of residues in the protein
	energy -- what the calculated energy of the protein is
	sequence -- requires use of the residue class, records references to 
		residues in order. Residues themselves contain their own interactive
		properties and the direction of the next residue relative to itself.
	proteinMap -- requires use of the coordinate class.  Is a secondary 
		construction to the sequence as it actually contains the coordinates
		that residues in the sequence would occupy given the relative 
		directions contained in sequence residue's attributes
	
	CONSTRUCTORS:
	Protein() -- the No-Argument constructor
	Protein(Residue firstResidue) -- Init Constructor.  this constructor is used when 
		a brand new protein is started.  The first residue is specified and 
		made the N-terminus, which is biologically the first terminus.
	Protein(Protein templateProtein) -- Copy Constructor.  makes a copy of the 
		protein in the parameters.
	
	METHODS:
	GETTERS
	SETTERS
	updateEnergy(): void -- recalculates the energy for the current sequence of
		residues, called when getEnergy is called to ensure current energy given.
	addResidue(Direction direction, ResidueType residue): void -- adds a new residue to the sequence
		and proteinMap. Needs to have the direction (part of the Direction enumerated data Type
		of the residue specified and should be determined by the decision making part of protein_Folder,
		the Colony class. The type of residue must also be given and must be
		from the enumerated data type ResidueType.
	removeResidue() : void -- removes the last residue in the chain, the C-terminus
	getResidue(int index): Residue -- The index of the residue, from the N-terminus, is given
		and returns a reference to the Residue object at that position.
	updateMap(): void -- Clears the current proteinMap and updates it according to current sequence
		by calling the actual mapping function mapProtein.
	mapProtein(): ArrayList<Coordinate> -- Reads the directions and types of residues from the sequence
		and associates each residue with a coordinate.  Returns a list of references to Coordinate
		objects in order of each sequence.  Each Coordinate object contains the axis coordinates and the
		residue type so that it can be printed easier.
	printCoordMap(): void -- simply prints the list of coordinates with the residue type.
	
	method prototypes
		public int getLength();
		public double getTotalEta();
		public double calcEta(Residue residue);
		public int getEnergy();
		public void updateEnergy();
		public void updateMap();
		
		public void addResidue(Direction direction); //adds a residue to the end of the chain in the specified direction
		public void changeResidue(Residue residue, Direction newDirection); //changes the direction of a specified residue
		public boolean checkAttrition();
		public void fixAttrition();  
		
		//output methods
		public Direction getResidue(int index);  //tells you the direction of residue index
		public void makeSequence();  //creates a text file of the sequence of residues
		public void makeProtein();  //create a text file of the protein
		public      mapProtein();  //returns a List(?) of the coordinates of each residue
					//the idea is that this could be used in checking for attrition, printing the protein out, and calculating the eta value
		
	
	Version 1.0 beta 
		- energy currently only takes into account the number of hydrogen bonds
		H-H between hydrophilic residues
		- Currently only 2D proteins are supported. Directions are L-left, R - right, and S - straight, 
		but extensibility is possible
		-no-arg constructor is currently not used
		-current ResidueTypes are H -hydrophobic,  P - polar and NULL
  ------------------------------------------------------------------------------
 */


package protein_folder;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.lang.Math;


import java.awt.*;
import java.awt.geom.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;




import javax.imageio.ImageIO;

import protein_folder.Param.Direction;
import protein_folder.Param.TerminalType;
import protein_folder.Param.UnitVector;
import protein_folder.Param.ResidueType;

public class Protein implements Comparable<Protein> {

		
//==================ATTRIBUTES======================================================
		int length; //number of residues in the protein
		int energy; //the number of H-H bonds in the protein
		double eta; //thermodynamic favorability of this protein
		
		LinkedList<Residue> sequence = new LinkedList<Residue>();  //the sequence of residues a protein currently has
		ArrayList<Coordinate> proteinMap = new ArrayList<Coordinate>();  //a map of coordinates for how the protein is currently folded
//==================================================================================	
		
//////__________Constructors_________________________________________________________
		
//----No-Arg Constructor-----------------------------------------------------------
		
		public Protein(){
			sequence.add(new Residue());
			proteinMap.add(new Coordinate());
			length = 0;
			energy = 0;
			eta = 0;
		}
//----End No-Arg Constructor-------------------------------------------------------
		
//----Init Constructor.-----------------------------------------------------------
		//constructor makes the N terminus, transcription occurs N to C termini
		public Protein(ResidueType inType){
			sequence.add(new Residue(inType, Param.N_TERMINUS));
			length = 1;
			eta = 1;
		}
//----End Init Constructor-------------------------------------------------------
		
//----Copy Constructor-----------------------------------------------------------

		public Protein(Protein inProt){
			
			//first copy basic attributes
			this.setLength(inProt.getLength());
			this.setEnergy(inProt.getEnergy());
			this.setEta(inProt.getEta());
			//then copy the arraylist attributes
			proteinMap.addAll(inProt.getMap());
			sequence.addAll(inProt.getSequence());
		}
//----End Copy Constructor-------------------------------------------------------
		
//___________________________________________________________________________________
		
//////++++++++GETTERS & SETTERS+++++++++++++++++++++++++++++++++++++++++++++++++++++
		public int getLength(){
			return length;
		}
		public void setLength(int inLength){
			this.length = inLength;
		}
		
		public int getEnergy(){
			updateEnergy();
			return energy;
		}
		public void setEnergy(int inEnergy){
			this.energy = inEnergy;
		}

		public double getEta() {
			return eta;
		}
		public void setEta(double eta) {
			this.eta = eta;
		}

		//get map function
		public ArrayList<Coordinate> getMap(){
			updateMap();
			return proteinMap;
		}
		//get residue sequence function
		public LinkedList<Residue> getSequence(){
			return sequence;
		}
		

		//get residue function
		public Residue getResidue(int index){
			return sequence.get(index);
		}
//////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	

		
//==================METHODS======================================================
		
		
//----makeImage-----------------------------------------------------------
static public BufferedImage makeImage(Protein prot, PNG_Footer footer) throws Exception{
	
	int width = 0, height = 0;
	int circleDiameter = Param.RESIDUE_CIRCLE_DIAMETER;
	int latticeSize = Param.LATTICE_SIZE;
	int footerSize = Param.FOOTER_SIZE;
	//the working map
	ArrayList<Coordinate> currProtMap = new ArrayList<Coordinate>();
	currProtMap.addAll(prot.getMap());
	
	//set the picture width and height for this conformation
	//find the min and max positions
	int maxX = 0, minX = 0, maxY = 0, minY = 0;			
	for(int i=0; i < currProtMap.size(); i++){
		//maxxes
		if(currProtMap.get(i).getX() > maxX)
			maxX = currProtMap.get(i).getX();
		if(currProtMap.get(i).getY() > maxY)
			maxY = currProtMap.get(i).getY();
		//mins
		if(currProtMap.get(i).getX() < minX)
			minX = currProtMap.get(i).getX();
		if(currProtMap.get(i).getY() < minY)
			minY = currProtMap.get(i).getY();
		
	}//end find max, min coordinates
	
	//now set the size of the image to be outputted
	width = (maxX-minX)*latticeSize + circleDiameter + Param.MARGINS*2;
	height = (maxY-minY)*latticeSize + circleDiameter + footerSize + Param.MARGINS*2 ;
	//a footer size is added to the bottom to display info about the conformation
	
	
	//the image
	BufferedImage bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
	//the graphics object
	Graphics2D ig2 = bi.createGraphics();
	//circle stroke
	BasicStroke circleOutline = new BasicStroke(Param.RESIDUE_CIRCLE_STROKE_WIDTH);
	//connector width set in a new stroke object
	BasicStroke connectorLine = new BasicStroke(Param.CONNECTOR_WIDTH);
	
	//make a white background
	ig2.setColor(Param.BACKGROUND_COLOR);
	ig2.fill(new Rectangle2D.Double(0,0,width,height));
	
	//draw the shapes for the map
	for(int i=0; i < currProtMap.size(); i++){
		//extract the coordinates for the current sphere
		int currX = currProtMap.get(i).getX(), currY = currProtMap.get(i).getY();
		//transform to the rendering space
		currX = (currX - minX) * latticeSize + Param.MARGINS;
		currY = (currY - minY) * latticeSize + Param.MARGINS;
		
		//draw the connector lines
		ig2.setStroke(connectorLine);
		ig2.setColor(Param.CONNECTOR_COLOR);
		
		if(i>0){
			//extract the coordinates for the past sphere
			int lastX = currProtMap.get(i-1).getX(), lastY = currProtMap.get(i-1).getY();
			//transform to the rendering space
			lastX = (lastX - minX) * latticeSize + Param.MARGINS;
			lastY = (lastY - minY) * latticeSize + Param.MARGINS;
			
			ig2.draw(new Line2D.Double(lastX + circleDiameter/2, lastY + circleDiameter/2,
					currX + circleDiameter/2, currY + circleDiameter/2) );
			//REDRAW the old circle so it is on top
			//draw the appropriate sphere
			ig2.setStroke(circleOutline);
			
			if(currProtMap.get(i-1).getType() == Param.POLAR){
				ig2.setColor(Param.STROKE_COLOR);
				ig2.draw(new Ellipse2D.Double(lastX, lastY, circleDiameter, circleDiameter) );
				ig2.setColor(Param.POLAR_COLOR);
				ig2.fill(new Ellipse2D.Double(lastX, lastY, circleDiameter, circleDiameter) );
			}
			else if(currProtMap.get(i-1).getType() == Param.HYDROPHOBIC){
				ig2.setColor(Param.STROKE_COLOR);
				ig2.draw(new Ellipse2D.Double(lastX, lastY, circleDiameter, circleDiameter) );
				ig2.setColor(Param.HYDROPHOBIC_COLOR);
				ig2.fill(new Ellipse2D.Double(lastX, lastY, circleDiameter, circleDiameter) );
			}
			
		}//end guard for connector drawing need to be on node 1 or more
		
		//draw the appropriate sphere
		ig2.setStroke(circleOutline);
		
		if(currProtMap.get(i).getType() == Param.POLAR){
			ig2.setColor(Param.STROKE_COLOR);
			ig2.draw(new Ellipse2D.Double(currX, currY, circleDiameter, circleDiameter) );
			ig2.setColor(Param.POLAR_COLOR);
			ig2.fill(new Ellipse2D.Double(currX, currY, circleDiameter, circleDiameter) );
		}
		else if(currProtMap.get(i).getType() == Param.HYDROPHOBIC){
			ig2.setColor(Param.STROKE_COLOR);
			ig2.draw(new Ellipse2D.Double(currX, currY, circleDiameter, circleDiameter) );
			ig2.setColor(Param.HYDROPHOBIC_COLOR);
			ig2.fill(new Ellipse2D.Double(currX, currY, circleDiameter, circleDiameter) );
		}
			
	}//end print spheres for loop
	
	//set the footer text
	Font font = new Font("TimesRoman", Font.BOLD, Param.FONT_SIZE);
	ig2.setFont(font);
	FontMetrics fontMetrics = ig2.getFontMetrics();
	//int stringWidth = fontMetrics.stringWidth(message);
	int stringHeight = fontMetrics.getAscent();
	ig2.setPaint(Param.FONT_COLOR);
	ig2.drawString(footer.getSolutionStr(), Param.MARGINS, height-footerSize + stringHeight );
	ig2.drawString(footer.getEnergyStr(), Param.MARGINS, height-footerSize + stringHeight*2 );
	ig2.drawString(footer.getSequenceInfoStr(), Param.MARGINS, height-footerSize + stringHeight*3 );
	ig2.drawString(footer.getColInfoStr(), Param.MARGINS, height-footerSize + stringHeight*4 );

	return bi;
	
	
}//end makeImage
	
	
//----End makeImage-------------------------------------------------------
		
//----updateEnergy-----------------------------------------------------------
		
		public void updateEnergy(){
			ListIterator<Coordinate> iter1 = proteinMap.listIterator();
			ListIterator<Coordinate> iter2 = proteinMap.listIterator();
			Coordinate res1 = new Coordinate();
			Coordinate res2 = new Coordinate();
			int seqPos1 = 0; //used to make sure adjacent residues on the chain dont get counted as an interaction
			int seqPos2 = 0;// 
			int H_Count = 0; // number of H-H bonds found in the protein
			
				//find H-H pairs
				while( iter1.hasNext()){
					//get the next residue for iter1
					res1 = iter1.next();
					
					//if this residue is polar we can look for a paired residue
					if(res1.getType() == Param.POLAR){
						//now compare the first residue to each other polar residue
						//make an H-H count if the two are both polar and non-sequentially adjacent 
						
						while(iter2.hasNext()){
							//get next residue for iter2
							res2 = iter2.next();
							
							
							//if non-sequentially adjacent, identical, and polar continue to test for adjacency
							if(res2.getType() == Param.POLAR &&
								seqPos2 != seqPos1 && 
								seqPos2 != seqPos1+1 && 
								seqPos2 != seqPos1-1){
								//if they are adjacent by one unit but not diagonally by using the distance formula
								if(Math.sqrt( Math.pow((double)res2.getX()-(double)res1.getX(), 2.0) + Math.pow((double)res2.getY()-(double)res1.getY(), 2.0) ) == 1 )
								{  
									H_Count++;  // add +1 to H-H count
								}//end adjacent if
								
							}//end non-sequential adjacency if
							seqPos2++;
						}//end end iter2 loop
						//reset position and iterator
						seqPos2 = 0;
						iter2 = proteinMap.listIterator();
						
					}//end if for iter1 residue
					seqPos1++;
				}//end iter1 loop
				
			
			//for each H-H bond subtract 1 from the total energy
				//divide H_Count by 2 because for each bond it will counted twice for each residue involved
			energy = 0 - H_Count/2;
			
		}//end update energy method
//----End updateEnergy-------------------------------------------------------
		
//----Function-----------------------------------------------------------	
		
		public void addResidue(Direction dir, ResidueType inType){
			//the direction relative to the existing chain is given by the residue previous it so 
			//this must be set 
			sequence.getLast().setNextResDirection(dir);
			
			//then we must remove its status as the C terminus unless it is the N terminus
			if(sequence.getLast().getTerminus() != Param.N_TERMINUS){
				sequence.getLast().setTerminus(Param.NOT_TERMINUS);
			}
			
			//add a new residue making it the C terminus, transcription is from N to C termini
			//with the null direction E for end
			sequence.add(new Residue(Param.END, inType, Param.C_TERMINUS));
			//increase the length attribute
			length++;
			this.updateMap();
			this.updateEnergy();
			
		}//end add residue
//----End addResidue-------------------------------------------------------
		
//----removeResidue-----------------------------------------------------------
		
		//the remove residue function removes the end residue the C-terminus
		//returns true if succesful false if not
		public boolean removeResidue(){
			//check to see it has >1 residue to remove
			if(sequence.size() > 1){
				sequence.removeLast();
				sequence.getLast().setNextResDirection(Param.END);
				sequence.getLast().setTerminus(Param.C_TERMINUS);
				this.length--;
				this.updateEnergy();
				return true;
			}else{
				System.out.println("Error: There is only one residue left");
				return false;
			}
			
			
		}//end the remove residue function
//----End removeResidue-------------------------------------------------------
		
//----updateMap---------------------------------------------------------------
		
		//updates the Map (coordinates of residues in lattice space) attribute by calling mapProtein
		public void updateMap(){
			proteinMap.clear();
			proteinMap.addAll(this.mapProtein());
		}
//----End updateMap-------------------------------------------------------
		
//----mapProtein-----------------------------------------------------------

		public ArrayList<Coordinate> mapProtein(){
			//initialize the straight (S) direction
				UnitVector str = Param.I_VECTOR;
			//Create an ArrayList<Coordinate> for storing the sequence of coordinates
				ArrayList<Coordinate> resMap = new ArrayList<Coordinate>();
			//define an iterator for the sequence List
				ListIterator<Residue> iterSeq = sequence.listIterator();
			//set the current residue to the first residue
				Residue currRes = sequence.getFirst();
			//set a variable to hold the current residue type
				ResidueType currType = currRes.getType();
			//create another temporary coordinate object to calculate the next coordinate from
				Coordinate oldCoord = new Coordinate();
				//is set to the initial coordinate already (0,0)
				
			
			//the first residue also hold the direction of the next residue so we get the second coord
			//use a do-while structure to implement this
			while(iterSeq.hasNext()){  //end while loop when last residue is reached in sequence list){
				
				//now update the unit vector and increment the coordinate for the next residue
				//move iterator to next residue
				currRes = iterSeq.next();
				//update current residue type
				currType = currRes.getType();
				
				//make a new coordinate object to put into the ArrayList 
				Coordinate newCoord = new Coordinate();
				//set the coordinate attributes
				newCoord.setX(oldCoord.getX());
				newCoord.setY(oldCoord.getY());
				newCoord.setType(currType);
				//now set the newCoord to the ArrayList
				resMap.add(newCoord);
				
				
				//set the unit vector direction the next residue will be mapped
				switch(currRes.getNextResDirection()){
					//next residue is straight ahead
					case S:  switch(str){ case i : 	//str is same
													break;
										  case I : //str is same
											  		break;
										  case j : //str is same
											  		break;
										  case J : //str is same
											  		break;
										  default: System.out.println ("invalid unit vector used");
										  		break;
										}
					         break;
					case L:  switch(str){ case i : str = Param.j_VECTOR; //y direction is straight
													break;
										  case I : str = Param.J_VECTOR; //-y direction is straight
													break;
										  case j : str = Param.I_VECTOR; //-x direction is straight
													break;
										  case J : str = Param.i_VECTOR; //x direction is straight
													break;
										  default: System.out.println ("invalid unit vector used");
										  			break;
										}
							break;
					case R:  switch(str){ case i : str = Param.J_VECTOR; //-y direction is straight
												break;
										  case I : str = Param.j_VECTOR;  //y direction is straight
												break;
										  case j : str = Param.i_VECTOR;  //x direction is straight
												break;
										  case J : str = Param.I_VECTOR;  //-x direction is straight
												break;
										  default: System.out.println ("invalid unit vector used");
										  			break;
										}
							break;
					case E: 	//System.out.println ("end of chain");
							break;
					default:  System.out.println ("error");
							break;
						     
			}
				
			//determine which coordinate must be incremented
			  //for the direction given above
				 //relative to the straight direction 				
					switch(str){  case i : oldCoord.incX();
									break;
								  case I : oldCoord.decX();
									break;
								  case j : oldCoord.incY();
									break;
								  case J : oldCoord.decY();
								   //str is same
									break;
								  default: System.out.println ("invalid unit vector used");
								  	break;
					}//end the coordinate setting
											
			}
			
			return resMap;
		}//end the mapProtein method
//----End mapProtein-------------------------------------------------------
		
//----printCoordMap-----------------------------------------------------------

		public void printCoordMap(){
			
			ArrayList<Coordinate> activeMap = new ArrayList<Coordinate>();
			activeMap.addAll(proteinMap);
			
			//print out map
			ListIterator<Coordinate> iterCoord = activeMap.listIterator();
			while(iterCoord.hasNext()){
				
				Coordinate currentCoord = iterCoord.next();
				System.out.println(currentCoord.stringCoord());
			}//end output loop
		}//end printCoordmap
		
//----End printCoordMap-------------------------------------------------------
		
//----Function compareTo-----------------------------------------------------------

	@Override
	public int compareTo(Protein compareProt) {
		
		int compareEnergy = ((Protein) compareProt).getEnergy(); 
		 
		//ascending order
		return this.energy - compareEnergy;
 
	}
		
//----End compareTo-------------------------------------------------------
	

		
//==================================================================================
		
}//end protein class
