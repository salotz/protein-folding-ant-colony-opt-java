package protein_folder;

import java.awt.image.BufferedImage;
import java.util.ArrayList;

public class ProteinBuild {

	
	private ArrayList<Protein> buildSteps;
	private ArrayList<BufferedImage> buildStepsImages;
	private String antType;
	private int wave;
	private int antNum;
	

	public ProteinBuild(){
		
		buildSteps = new ArrayList<Protein>();
		buildStepsImages = new ArrayList<BufferedImage>();
		wave=0;
		antNum = 0;
		
		
	}//end no-arg constructor
	
	public ArrayList<Protein> getBuildSteps() {
		return buildSteps;
	}
	public void setBuildSteps(ArrayList<Protein> buildSteps) {
		this.buildSteps = buildSteps;
	}
	public ArrayList<BufferedImage> getBuildStepsImages() {
		return buildStepsImages;
	}
	public void setBuildStepsImages(ArrayList<BufferedImage> buildStepsImages) {
		this.buildStepsImages = buildStepsImages;
	}
	public String getAntType() {
		return antType;
	}
	public void setAntType(String antType) {
		this.antType = antType;
	}
	public int getWave() {
		return wave;
	}
	public void setWave(int wave) {
		this.wave = wave;
	}
	public int getAntNum() {
		return antNum;
	}
	public void setAntNum(int antNum) {
		this.antNum = antNum;
	}
	
}//end class
