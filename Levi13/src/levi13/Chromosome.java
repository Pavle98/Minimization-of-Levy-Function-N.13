package levi13;

public class Chromosome implements Comparable<Chromosome> {
	public double x;
	public double y;

	
	public Double fitnessFunction;


	@Override
	public int compareTo(Chromosome o) {
		return this.fitnessFunction.compareTo(o.fitnessFunction);
	}
}