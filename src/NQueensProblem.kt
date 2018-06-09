import kotlin.math.abs

fun main(args: Array<String>) {
    val queens = NQueensProblem()

    queens.makePopulation()

    println("=========================================================")
    println("INITIAL GENERATION:")
    printPopulation(queens)
    queens.evalPopulation(queens.population)
    printFitness(queens)
    queens.getSummary(0)
    queens.makeGenerations()
}

private fun printPopulation(kp: NQueensProblem) {
    println("\nPopulation:")
    kp.population.forEachIndexed { index, elem -> println("${index + 1} - $elem") }
}

private fun printFitness(kp: NQueensProblem) {
    println("\nFitness:")
    kp.fitness.forEachIndexed { index, elem -> println("${index + 1} - $elem (${1000/elem} conflicts)") }
}

class NQueensProblem {
    /**
     * Setup
     */
    private var n = 10
    private var populationSize: Int = 20
    private var maxGenerations: Int = 600
    private val crossoverProbability: Double = 0.5
    private val mutationProbability: Double = 0.3

    /**
     * Init variables
     */
    val population: MutableList<String> = mutableListOf()
    private var totalGenerationFitness: Double = 0.0
    var fitness: MutableList<Double> = mutableListOf()
    private val bestGenerationSolution: MutableList<String> = mutableListOf()
    private val bestGenerationFitness: MutableList<Double> = mutableListOf()
    private val averageGenerationFitness: MutableList<Double> = mutableListOf()
    private var mutationCount: Int = 0
    private var crossoverCount: Int = 0
    private var cloneCount: Int = 0
    private var generationCounter: Int = 1
    private val breedPopulation: MutableList<String> = mutableListOf()

    /**
     * Filling population with random genes
     */
    fun makePopulation() {
        for (i in 0 until populationSize) {
            population.add(generateGene())
        }
    }

    /**
     * Generates gene - random queens placement on the cheesboard
     */
    private fun generateGene(): String {
        var gene = ""
        for (i in 0 until n) {
            gene += (1..n + 1).random().toString()
        }
        return gene
    }

    /**
     * Evaluates population fitness
     */
    fun evalPopulation(population: MutableList<String>) {
        totalGenerationFitness = 0.0
        for (i in 0 until populationSize) {
            val tmpFitness = evalGene(population[i])
            fitness.add(tmpFitness)
            totalGenerationFitness += tmpFitness
        }
    }

    /**
     * Evaluates a single gene's fitness
     */
    private fun evalGene(gene: String): Double {
        val chessboard = hashMapOf<Int, Char>()
        var fitnessValue = 0.0
        for (i in 0 until n) {
            chessboard[i + 1] = gene[i]
        }

        val matchedPairs: MutableList<Pair<Pair<Int, Char>, Pair<Int, Char>>> = mutableListOf()
        for (i in 0 until n) {
            val currentX = i + 1
            val currentY = chessboard[i + 1]
            for (j in 0 until n) {
                val tmpX = j + 1
                val tmpY = chessboard[j + 1]
                if (isConflictingAcross(tmpX, currentX, tmpY, currentY, matchedPairs)) {
                    matchedPairs.add(Pair(Pair(currentX, currentY!!), Pair(tmpX, tmpY!!)))
                    fitnessValue += 1
                }
                if (abs(currentX - tmpX) == abs(currentY.toString().toInt() - tmpY.toString().toInt()) &&
                    (currentX != tmpX && currentY != tmpY) &&
                    !matchedBefore(matchedPairs, tmpX, tmpY, currentX, currentY)
                ) {
                    matchedPairs.add(Pair(Pair(currentX, currentY!!), Pair(tmpX, tmpY!!)))
                    fitnessValue += 1
                }
            }
        }
        return 1 / fitnessValue * 1000
    }

    private fun isConflictingAcross(
        tmpX: Int,
        currentX: Int,
        tmpY: Char?,
        currentY: Char?,
        matchedPairs: MutableList<Pair<Pair<Int, Char>, Pair<Int, Char>>>
    ): Boolean {
        return tmpX != currentX && tmpY == currentY && !matchedBefore(matchedPairs, tmpX, tmpY, currentX, currentY)
    }

    private fun matchedBefore(
        matchedPairs: MutableList<Pair<Pair<Int, Char>, Pair<Int, Char>>>,
        tmpX: Int,
        tmpY: Char?,
        currentX: Int,
        currentY: Char?
    ): Boolean {
        return matchedPairs.contains(
            Pair(
                Pair(tmpX, tmpY),
                Pair(currentX, currentY)
            )
        )
    }

    private val bestSolution: Int
        get() {
            var bestPosition = 0
            var currentFitness = 0.0
            var bestFitness = 0.0
            for (i in 0 until populationSize) {
                currentFitness = evalGene(population[i])
                if (currentFitness > bestFitness) {
                    bestFitness = currentFitness
                    bestPosition = i
                }
            }
            return bestPosition
        }

    private val averageFitness: Double
        get() {
            var totalFitness = 0.0
            var averageFitness = 0.0
            for (i in 0 until populationSize) {
                totalFitness += fitness[i]
            }
            averageFitness = totalFitness / populationSize
            return averageFitness
        }

    fun getSummary(generationIndex: Int) {
        bestGenerationSolution.add(population[bestSolution])
        bestGenerationFitness.add(evalGene(population[bestSolution]))
        averageGenerationFitness.add(averageFitness)
        println("\n-----------------------------------")
        println("| Best solution: ${bestGenerationSolution[generationIndex]}")
        println("| Best fitness: ${bestGenerationFitness[generationIndex]} (${1000/bestGenerationFitness[generationIndex]} conflicts)")
        println("| Average fitness: ${averageGenerationFitness[generationIndex]}")
        println("-----------------------------------")
        println("| Crossover:  $crossoverCount times")
        println("| Cloning:  $cloneCount times")
        println("| Mutation:  $mutationCount times")
        println("-----------------------------------")
    }

    private fun tournamentSelection(): String {
        val randPickOne = (0..populationSize).random()
        val randPickTwo = (0..populationSize).random()
        val randPickThree = (0..populationSize).random()

        val topGene = mutableListOf(
            Pair(population[randPickOne], fitness[randPickOne]),
            Pair(population[randPickTwo], fitness[randPickTwo]),
            Pair(population[randPickThree], fitness[randPickThree])
        )
            .sortedByDescending { it.second }
            .first()

        return topGene.first
    }

    fun makeGenerations() {

        for (i in 1 until maxGenerations) {
            if (checkForStopCriteria(i)) break
            resetCounters()

            while (breedPopulation.size < populationSize) {
                // if 2 genes wont fit into new population just copy best solution from previous generation
                if (populationSize - breedPopulation.size == 1)
                    breedPopulation.add(bestGenerationSolution[generationCounter - 1])
                crossoverGenes(tournamentSelection(), tournamentSelection())
                mutateGene()
            }

            // Clear fitness values of previous generation
            fitness.clear()

            // Evaluate fitness of breed population members
            evalPopulation(breedPopulation)

            // Copy breed_population to population
            for (k in 0 until populationSize) {
                population[k] = breedPopulation[k]
            }

            println("=========================================================")
            println("\nGENERATION ${(i + 1)}")
            printPopulation(this)
            printFitness(this)
            breedPopulation.clear()
            getSummary(i)
        }
    }

    private fun resetCounters() {
        crossoverCount = 0
        cloneCount = 0
        mutationCount = 0
    }

    private fun checkForStopCriteria(i: Int): Boolean {
        if (bestGenerationFitness[i-1].isInfinite()) {
            println("Found solution")
            return true
        }

        if (maxGenerations > 4 && i > 4) {

            // Previous 3 generational average fitness values
            val a = averageGenerationFitness[i - 1]
            val b = averageGenerationFitness[i - 2]
            val c = averageGenerationFitness[i - 3]

            if (a == b && b == c) {
                println("\nStop criteria!")
                return true
            }
        }
        return false
    }

    /**
     * Crossover
     */
    private fun crossoverGenes(geneOne: String, geneTwo: String) {

        val crossoverRand = Math.random()
        if (crossoverRand <= crossoverProbability) {
            crossoverCount += 1
            val crossingPoint = (0..n).random()

            // Crossing genes at randomly chosen spot
            val newGeneOne = geneOne.substring(0, crossingPoint) + geneTwo.substring(crossingPoint)
            val newGeneTwo = geneTwo.substring(0, crossingPoint) + geneOne.substring(crossingPoint)

            // Add new genes to breed_population
            breedPopulation.add(newGeneOne)
            breedPopulation.add(newGeneTwo)
        } else {
            cloneCount += 1
            breedPopulation.add(geneOne)
            breedPopulation.add(geneTwo)
        }
    }

    /**
     * Mutation
     */
    private fun mutateGene() {
        val mutationRand = Math.random()
        if (mutationRand <= mutationProbability && breedPopulation.size >= 1) {
            mutationCount += 1
            val mutatedGene: String = breedPopulation[(0..breedPopulation.size).random()]
            val mutationPoint = (0..n).random()
            val newGene = mutatedGene.replaceRange(mutationPoint, mutationPoint + 1, (1..n + 1).random().toString())

            breedPopulation[breedPopulation.indexOf(mutatedGene)] = newGene
        }
    }
}
