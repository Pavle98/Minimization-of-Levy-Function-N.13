package levi13;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.Math;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
// Zadatak je Minimizacija Levijeve funkcije broj 13 na intervalu -10<= x,y <= 10, ocekivani minimum f(1,1) = 0
// Tip genetskog algoritma = Kontinualni genetski algoritam
//Metoda ukrstanja = BLX-α ukrstanje za kontinualni genetski algoritam
//Metoda mutacije = Tackasta uniformna mutacija za kontinualni genetski algoritam ( bira se slucajna vrednost iz domena gena po uniformnoj raspodeli)
public class Main {
	static int repetitions = 3;	//maximum number of experiment repetitions
	public static int maxIterations = 500; //maximum number of generations in one repetition
	static double mutationRate = 0.2; // chance of mutation, so.... 20% here...
	public static double alpha = (float) 0.5; //alpha value for blx
	static int populationNumber = 100; // number of chroms in one population
	static int[] scope = { -10, 10 }; //  (given interval)
	static int tournamentPopulationNumber = 4; //numbe of chroms that participate in tournamate in selection 

	public static void main(String[] args) {
		myalgorithm();
	}
	//Deffinition of LEVY FUNCTION N. 13
	public static double fitness_function(Chromosome chroms) {
		double x = chroms.x;
		double y = chroms.y;
		return Math.pow(Math.sin(3 * Math.PI * x), 2)
				+ Math.pow(x - 1, 2) * (1 + Math.pow(Math.sin(3 * Math.PI * y), 2))
				+ Math.pow(y - 1, 2) * (1 + Math.pow(Math.sin(2 * Math.PI * y), 2));
	}
	//Uniform mutation that uses a random value from the scope
	public static Chromosome mutate(Chromosome chroms, double mut_rate, int[] scope) {
		double rand = (Math.random());
		if (rand < mut_rate) {
			chroms.x = scope[0] + (Math.random() * scope[1] * 2);
			chroms.y = scope[0] + (Math.random() * scope[1] * 2);
		}
		return chroms;
	}

	//
	public static Chromosome selection(ArrayList<Chromosome> population, int size) {
		ArrayList<Chromosome> tournamentList = new ArrayList<Chromosome>();
		Random rand = new Random();
		int randChrom = rand.nextInt(population.size());
		double selectedValue = 0;
		Chromosome bestChrom = null;
		double bestFunctionValue = -100; /// ovo vrv nije dobro

		while (tournamentList.size() < population.size()) {
			tournamentList.add(population.get(randChrom));
			//System.out.println(randChrom + " rando chrom ");
		}

		for (Chromosome selectedChrom : tournamentList) {
			selectedValue = fitness_function(selectedChrom);
		
			if (bestChrom == null || selectedValue > bestFunctionValue) {
				bestFunctionValue = selectedValue;
				bestChrom = selectedChrom;
			}

		}
		return bestChrom;
	}
	//BLX-alpha crossover
	public static ArrayList<Chromosome> crossover(Chromosome firstChrom, Chromosome secondChrom) {
		ArrayList<Chromosome> newChroms = new ArrayList<Chromosome>();
		Chromosome firstNewChrom = new Chromosome();
		Chromosome SecondNewChrom = new Chromosome();
		newChroms.add(firstNewChrom);
		newChroms.add(SecondNewChrom);
		Random rand = new Random();
		double d;
		for (int i = 0; i < 2; i++) {
			if (i == 0) {
				d = Math.abs(firstChrom.x - secondChrom.x);
				double min = Math.min(firstChrom.x, secondChrom.x) - (alpha * d);
				double max = Math.max(firstChrom.x, secondChrom.x) + (alpha * d);
				double u = rand.nextDouble() * (max - min) + min;
				firstNewChrom.x = u;

				min = Math.min(firstChrom.x, secondChrom.x) - (alpha * d);
				max = Math.max(firstChrom.x, secondChrom.x) + (alpha * d);
				u = rand.nextDouble() * (max - min) + min;
				SecondNewChrom.x = u;
			} else {
				d = Math.abs(firstChrom.y - secondChrom.y);
				double min = Math.min(firstChrom.y, secondChrom.y) - (alpha * d);
				double max = Math.max(firstChrom.y, secondChrom.y) + (alpha * d);
				double u = rand.nextDouble() * (max - min) + min;
				firstNewChrom.y = u;

				min = Math.min(firstChrom.y, secondChrom.y) - (alpha * d);
				max = Math.max(firstChrom.y, secondChrom.y) + (alpha * d);
				u = rand.nextDouble() * (max - min) + min;
				SecondNewChrom.y = u;
			}

		}
		return newChroms;

	}

	public static void myalgorithm() {
		 BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter("Chroms.txt"));
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		System.out.println("Alg starting");
		int newPopulationSize = populationNumber;
		Chromosome bestChromosome;
		ArrayList<Chromosome> initialPopulation = new ArrayList<Chromosome>();
		ArrayList<Chromosome> newPopulation = new ArrayList<Chromosome>();
		double bestFitness = -100;
		int iter = 0;
		for (int i = 0; i < repetitions; i++) {
			initialPopulation = new ArrayList<Chromosome>();
			iter = 0;
			System.out.println("repetition " + i);
			bestChromosome = new Chromosome();
			bestChromosome.fitnessFunction = 1000.0;
			bestFitness = 100;

			for (int popi = 0; popi < populationNumber; popi++) {
				Chromosome andAnotherOne = new Chromosome();
				andAnotherOne.x = scope[0] + (Math.random() * scope[1] * 2);
				andAnotherOne.y = scope[0] + (Math.random() * scope[1] * 2);
				
				initialPopulation.add(andAnotherOne);
			}
			

			while (bestChromosome.fitnessFunction > 0.0000001 && iter < maxIterations && (bestChromosome.x != 1.0 && bestChromosome.y != 1.0 )) {

				newPopulation.removeAll(newPopulation);
					for(int ra = 0; ra <populationNumber; ra++) {
						newPopulation.add(initialPopulation.get(ra));
					}
				while(newPopulation.size() < (populationNumber*2)) {
				
					Chromosome chr1 = selection(initialPopulation, tournamentPopulationNumber);
					Chromosome chr2 = selection(initialPopulation, tournamentPopulationNumber);
					ArrayList<Chromosome> chrsFromCrossover = crossover(chr1, chr2);
					Chromosome chr3 = chrsFromCrossover.get(0);
					Chromosome chr4 = chrsFromCrossover.get(1);
					mutate(chr3, mutationRate, scope);
					mutate(chr4, mutationRate, scope);
					newPopulation.add(chr3);
					newPopulation.add(chr4);
				}
				for(Chromosome newPopChrom : newPopulation) {
					newPopChrom.fitnessFunction = fitness_function(newPopChrom);
				}
				Collections.sort(newPopulation);
				initialPopulation.removeAll(newPopulation);
				for(int r = 0; r <populationNumber; r++) {
					initialPopulation.add(newPopulation.get(r));
				}
				for(int r = 0; r <200; r++) {
				}
				double bestFitValue = fitness_function(newPopulation.get(0));
				
				if (bestFitness > bestFitValue) {
					bestFitness = bestFitValue;
					bestChromosome = newPopulation.get(0);
				}
				double totalfitness = 0;
				for(int ik = 0; ik < populationNumber ; ik++) {
					totalfitness+= newPopulation.get(ik).fitnessFunction;
				}
				double avg = totalfitness / newPopulation.size();
				
				
				String data = "Generation" + iter+":"+ " Best chromosome fitness =  " +bestFitness  + " Average  ="  + avg + "\n";
				System.out.printf("Generation " + iter +":"+ " Best chromosome fitness = %.10f "  + " Average  ="  + avg + "\n", bestFitness);
			    try {
			        writer.write(data);
		        } catch (Exception e) {
		            e.printStackTrace();
		        }
				
				iter += 1;
			}
			
			
			System.out.println("Best chromosome [" +bestChromosome.x +"] ["+ bestChromosome.y + "]");
			try {
				writer.write("Best chromosome [" +bestChromosome.x +"] ["+ bestChromosome.y + "] \n");		
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		try {
			
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
