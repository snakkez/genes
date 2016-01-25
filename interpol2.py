import random
from math import sqrt


class GeneticInterpolator(object):
    def __init__(self, settings):
        self.x_start, self.x_stop = settings["range_x"]
        self.y_start, self.y_stop = settings["range_y"]
        self.range_len = len(range(self.x_start, self.x_stop))
        self.chromosome_length = settings["chromosome_length"]
        if self.chromosome_length == 8:
            self.calc_func = self.func8
        elif self.chromosome_length == 6:
            self.calc_func = self.func6
        elif self.chromosome_length == 4:
            self.calc_func = self.func4
        else:
            raise Exception

        self.target_fit = settings["target_fit"]

        self.init_pop_size = settings["init_pop_size"]
        self.pop_cap = settings["pop_cap"]
        self.mutation_factor = settings["mutation_factor"]
        self.survival_factor = settings["survival_factor"]
        self.reproduction_factor = settings["reproduction_factor"]
        self.offspring_limit = settings["offspring_limit"]
        self.generation = 0
        self.global_specimen_counter = 0
        self.population = []
        self.average_score = 0
        self.lowest_score = 0
        target_chromosome = [1,1,1,1]
        self.target_result = self.calc_func(target_chromosome)

        for i in range(self.init_pop_size):
            self.population.append(Specimen(specimen_num=self.get_specimen_counter(),
                                            chromosome_length=self.chromosome_length,
                                            y_start=self.y_start, y_stop=self.y_stop))
        self.score_and_sort()

    def __str__(self):
        retstr= ""
        retstr += "Current generation: {:n}, average score: {: .3f}, lowest score: {: .3f}.\n".format(
            self.generation, self.average_score, self.lowest_score)
        retstr += "Surv.: {: .3f}, Mut: {: .3f}, Repr.: {: .3f},  OL: {: .3f}.\n".format(
            self.survival_factor, self.mutation_factor, self.reproduction_factor, self.offspring_limit)
        retstr += "Current population: {}\n".format(len(self.population))
        retstr += "-"*100+"\n"
        for specimen in self.population:
            retstr += str(specimen)+"\n"
        return retstr


    def func8(self, ch):
        return [ch[7]*(x**7) + ch[6]*(x**6) + ch[5]*(x**5) + ch[4]*(x**4) + ch[3]*(x**3) + ch[2]*(x**2) + ch[1]*x + ch[0]
                for x in range(self.x_start, self.x_stop)]
    def func6(self, ch):
        return [ch[5]*(x**5) + ch[4]*(x**4) + ch[3]*(x**3) + ch[2]*(x**2) + ch[1]*x + ch[0]
                for x in range(self.x_start, self.x_stop)]

    def func4(self, ch):
        return [ch[3]*(x**3) + ch[2]*(x**2) + ch[1]*x + ch[0] for x in range(self.x_start, self.x_stop)]


    def get_specimen_counter(self):
        self.global_specimen_counter += 1
        return self.global_specimen_counter

    def eval_fitness(self, specimen):
        tmp = []
        for i in range(self.range_len):
            calc_res = self.calc_func(specimen.chromosome)
            tmp.append((calc_res[i]-self.target_result[i]) ** 2)
        derivation = sqrt(sum(tmp)/len(tmp))
        specimen.score = derivation

    def score_and_sort(self):
        for specimen in self.population:
            self.eval_fitness(specimen)

        self.average_score = sum([specimen.score for specimen in self.population]) / len(self.population)

        self.population = sorted(self.population, key=lambda k: k.score)
        self.lowest_score = self.population[0].score

    def short_info(self):
        retstr= ""
        retstr += "Current generation: {:n}, average score: {: .3f}, lowest score: {: .3f}.\n".format(
            self.generation, self.average_score, self.lowest_score)
        retstr += "Surv.: {: .3f}, Mut: {: .3f}, Repr.: {: .3f},  OL: {: .3f}.\n".format(
            self.survival_factor, self.mutation_factor, self.reproduction_factor, self.offspring_limit)
        retstr += "Current population: {}\n".format(len(self.population))
        retstr += "-"*100+"\n"
        for ln in self.population[:5]:
            retstr += str(ln)+"\n"
        retstr+= "*"*100+"\n"
        for ln in self.population[-5:]:
            retstr += str(ln)+"\n"
        return retstr

    def breed(self, specimen1, specimen2):
        mix_point = random.randint(0, self.chromosome_length)
        offspring_chromosome = specimen1.chromosome[mix_point:]+specimen2.chromosome[:mix_point]
        if random.random() <= self.reproduction_factor:
            return Specimen(specimen_num=self.get_specimen_counter(), chromosome_length=self.chromosome_length,
                            inp_chromosome=offspring_chromosome)
        else:
            return None

    def mutate(self, specimen):
        new_chromosome = []
        for g in specimen.chromosome:
            if random.random() <= self.mutation_factor:
                g = random.uniform(self.y_start, self.y_stop)
            new_chromosome.append(g)
        specimen.chromosome = new_chromosome
        return specimen

    def survive(self, specimen):
        survival_probability = -(1/(self.survival_factor*self.average_score))*specimen.score + 1
        if random.random() <= survival_probability:
            return specimen
        else:
            return None


    def evolve(self):
        self.generation += 1
        previous_generation_pop = self.population
        self.population = []
        #Breeding
        pairs = [previous_generation_pop[i:i+2] for i in range(0, len(previous_generation_pop), 2)]
        for pair in pairs:
            if len(pair) == 2:
                self.population += [pair[0]] if pair[0] is not None else []
                self.population += [pair[1]] if pair[1] is not None else []
                for x in range(1, random.randrange(1, self.offspring_limit)):
                    offspring = self.breed(pair[0], pair[1])
                    self.population += [offspring] if offspring is not None else []
        #mutating
        self.population = [self.mutate(specimen) for specimen in self.population]


        #rescoring
        if len(self.population) > 0:
            self.score_and_sort()
            #capping
            if len(self.population) > self.pop_cap:
                self.population = self.population[:self.pop_cap]
            #surviving
            self.population = [self.survive(specimen) for specimen in self.population]
            self.population = [specimen for specimen in self.population if specimen is not None]



class Specimen(object):
    def __init__(self, specimen_num, chromosome_length, y_start=None, y_stop=None, inp_chromosome=None):
        if inp_chromosome is None:
            self.chromosome = []
            for i in range(chromosome_length):
                self.chromosome.append(random.uniform(y_start, y_stop))
        else:
            self.chromosome = inp_chromosome

        self.specimen_num = specimen_num
        self.chromosome_length = chromosome_length
        self.score = None

    def __repr__(self):
        return {"chromosome": self.chromosome, "score": self.score}

    def __str__(self):
        if self.score is not None:
            str_template = "#{:04d}:    "+"{: .3f}    " * self.chromosome_length + "    Score: {: .3f}"
            return str_template.format(*([self.specimen_num]+self.chromosome[::-1]+[self.score]))
        else:
            str_template = "#{:04d}:    "+"{: .3f}    " * self.chromosome_length + "    Score: N/A"
            return str_template.format(*([self.specimen_num]+self.chromosome[::-1]))






gp_settings = {
    "range_x": (0, 8),
    "x_step": 0.1,
    "range_y": (-500, 500),
    "chromosome_length": 4,
    "target_fit": 0.1,
    "mutation_factor": 0.1,
    "survival_factor": 4,
    "reproduction_factor": 1,
    "offspring_limit": 7,
    "init_pop_size": 500,
    "pop_cap": 2000
}

gp = GeneticInterpolator(gp_settings)

print(gp)
for x in range(1000):
    gp.evolve()
    print(gp.short_info())