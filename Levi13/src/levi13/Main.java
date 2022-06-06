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
		geneticAlgorithm();
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

	//selecting best chromosome from tournament of (size) random chromosomes
	public static Chromosome selection(ArrayList<Chromosome> population, int size) {
		ArrayList<Chromosome> tournamentList = new ArrayList<>();
		Random rand = new Random();
		int randChrom = rand.nextInt(population.size());
		double selectedValue = 0;
		Chromosome bestChrom = null;
		double bestFunctionValue = -100;

		while (tournamentList.size() < size) {
			tournamentList.add(population.get(randChrom));
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
	public static ArrayList<Chromosome> crossover(Chromosome firstParent, Chromosome secondParent) {
		ArrayList<Chromosome> newChroms = new ArrayList<>();
		Chromosome firstChild = new Chromosome();
		Chromosome secondChild = new Chromosome();
		newChroms.add(firstChild);
		newChroms.add(secondChild);
		Random rand = new Random();
		double d;
		double min;
		double max;
		double u;
		//Setting first and second children's X  based on BLX-a crossover
		d = Math.abs(firstParent.x - secondParent.x);
		min= Math.min(firstParent.x, secondParent.x) - (alpha * d);
		max = Math.max(firstParent.x, secondParent.x) + (alpha * d);
		u = rand.nextDouble() * (max - min) + min;
		firstChild.x = u;

		u = rand.nextDouble() * (max - min) + min;

		secondChild.x = u;

		//Setting Y
		d = Math.abs(firstParent.y - secondParent.y);
		min = Math.min(firstParent.y, secondParent.y) - (alpha * d);
		max = Math.max(firstParent.y, secondParent.y) + (alpha * d);
		u = rand.nextDouble() * (max - min) + min;
		firstChild.y = u;

		min = Math.min(firstParent.y, secondParent.y) - (alpha * d);
		max = Math.max(firstParent.y, secondParent.y) + (alpha * d);
		u = rand.nextDouble() * (max - min) + min;
		secondChild.y = u;

		return newChroms;

	}
	//Just filling population with chromosomes with random values from the scope
	public static void initialization (ArrayList<Chromosome> initialPopulation){

		for (int popi = 0; popi < populationNumber; popi++) {
			Chromosome andAnotherOne = new Chromosome();
			andAnotherOne.x = scope[0] + (Math.random() * scope[1] * 2);
			andAnotherOne.y = scope[0] + (Math.random() * scope[1] * 2);

			initialPopulation.add(andAnotherOne);
		}

	}

	public static void geneticAlgorithm() {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter("Chroms.txt"));
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		System.out.println("Alg starting");
		Chromosome bestChromosome;
		ArrayList<Chromosome> initialPopulation;
		ArrayList<Chromosome> newPopulation = new ArrayList<>();
		double bestFitness = -100;
		int iter = 0;
		for (int i = 0; i < repetitions; i++) {
			//ALgorithm repetitions, instead of starting program multiple times, we have
			//this for loop in which we reset initial parameters for each repetition
			// such as initial poppulation, best chrom and fitness, after that, we fill initialPop with random chroms
			initialPopulation = new ArrayList<>();
			iter = 0;
			System.out.println("repetition " + i);
			bestChromosome = new Chromosome();
			bestChromosome.fitnessFunction = 1000.0;
			bestFitness = 100;

			initialization(initialPopulation);

			//We are going through generations until we are satisfied with fitness function (7th decimal)
			//Or until we have done enough (specified number) of generations
			while (bestChromosome.fitnessFunction > 0.0000001 && iter < maxIterations) {

				newPopulation.removeAll(newPopulation);
				for(int ra = 0; ra <populationNumber; ra++) {
					newPopulation.add(initialPopulation.get(ra));
				}
				//after adding 100 of old chromosomes (From previous generation)
				//we are now adding 100 new chromosomes that have been through selection, crossover and mutation
				//thats why we have this condition
				while(newPopulation.size() < (populationNumber*2)) {

					//First we choose 2 chars via selection (Tournament)
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
				//We are sorting old + new chromosomes (Offsprings) by fitness value
				//(200 in total), where we will only choose top 100 for new generation
				for(Chromosome newPopChrom : newPopulation) {
					newPopChrom.fitnessFunction = fitness_function(newPopChrom);
				}

				Collections.sort(newPopulation);
				initialPopulation.removeAll(initialPopulation);
				for(int r = 0; r <populationNumber; r++) {
					initialPopulation.add(newPopulation.get(r));
				}

				//Calculating total, best and avg fitnesses, and writing them into file
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
