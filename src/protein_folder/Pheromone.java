package protein_folder;

import java.util.ArrayList;
import java.util.Iterator;

public class Pheromone {

	private int minX, maxX, minY, maxY;
	private ArrayList<Coordinate> resMap;
	private ArrayList<ArrayList<Double>> pheroMatrix;
	
	//METHODS
	//no arg constructor creates a pheromone matrix of default size 200x200
	Pheromone(int sizeX, int sizeY){
		
		minX= -(sizeX);
		maxX= sizeX;
		minY = -(sizeY);
		maxY = sizeY;
		
		//sets all spots to tau = 1.0
		for(int i=0; i < maxY-minY+1; i++){
			ArrayList<Double> newCol = new ArrayList<Double>();
			pheroMatrix.add(newCol);
			for(int j=0; j < maxX-minX+1; j++){
				Double newValue = new Double(Param.UNITY);
				pheroMatrix.get(i).add(newValue);
			}//end loop through the columns
			
		}//end loop through rows
	}//end no-arg constructor
	
	
	//argumented constructor for a given map creates a first pheromone matrix
	Pheromone(ArrayList<Coordinate> inResMap){
		
		resMap.addAll(inResMap);
		
		//find the max and min, x and y coordinates
		for(int seqCount =0; seqCount < resMap.size(); seqCount++){
			//min X
			if(minX > (int)resMap.get(seqCount).getX())
				minX = resMap.get(seqCount).getX();
			//maxX
			if(maxX < resMap.get(seqCount).getX())
				maxX = resMap.get(seqCount).getX();
			//min Y
			if(minY > resMap.get(seqCount).getY())
				minY = resMap.get(seqCount).getY();
			//maxY
			if(maxY < resMap.get(seqCount).getY())
				maxY = resMap.get(seqCount).getY();
		}//end min-max settings
		
		///DEBUG
		System.out.println("minX: " + minX);
		System.out.println("maxX: " + maxX);
		System.out.println("minY: " + minY);
		System.out.println("maxY: " + maxY);
		
		//initialize the pheroMatrix to 1 (UNITY)
		for(int i=0; i < maxY-minY+1; i++){
			ArrayList<Double> newCol = new ArrayList<Double>();
			pheroMatrix.add(newCol);
			for(int j=0; j < maxX-minX+1; j++){
				Double newValue = new Double(Param.UNITY);
				pheroMatrix.get(i).add(newValue);
			}//end loop through the columns
			
		}//end loop through rows
		
		
		
	}//end initial pheroMatrix constructor
	
	
	
	//print out the pheromone matrix
	public void printPheroMatrix(){
	
		//check to see if the pheromone matrix has been initialized
		if(resMap.isEmpty()){
			System.out.println("This pheromone matrix has not been initialized with a protein map");
		}
		else{
			
			System.out.println("Pheromone Matrix: "  );
			//iterate through the matrix to print out the values
			Iterator<ArrayList<Double>> row = pheroMatrix.iterator();
			
			while(row.hasNext()){
				
				Iterator<Double> col = row.next().iterator();
				
				while(col.hasNext()){
					System.out.print(col.next() + "  ");
				}//end column loop
				
				System.out.println();
				
			}//end row loop
			
			System.out.println();
		}//end check for existence
	}//end print pheromone matrix method
	
	//GETTERS AND SETTERS

	public ArrayList<ArrayList<Double>> getPheroMatrix() {
		return pheroMatrix;
	}


	public void setTau(int i, int j, Double newTau ){
		pheroMatrix.get(i).set(j, newTau);
	}


	
}//end pheromone class
