import unittest
from statistical_potential import *

class TestStatisticalPotential(unittest.TestCase):
    def test_observed_distribution_basic(self):
        train_structures = [
            [
                ('ASP', 'CB', 1, (0, 0, 0)),
                ('LYS', 'CB', 10, (0, 0, 1))
            ],
        ]

        x, y_obs = observed_distribution('ASP', 'LYS', train_structures)
        self.assertEqual(x, list(range(CUTOFF)))

        output = [EPS]* CUTOFF
        output[1] += 1
        self.assertEqual(y_obs, output)

    def test_observed_distribution_multiple_residues(self):
        train_structures = [
            [
                ('ASP', 'CB', 1, (0, 0, 0)),
                ('LYS', 'CB', 10, (0, 0, 1)),
                ('ARG', 'CB', 10, (0, 10.5, 0))
            ],
        ]

        x, y_obs = observed_distribution('ASP', 'ARG', train_structures)
        self.assertEqual(x, list(range(CUTOFF)))

        output = [EPS]* CUTOFF
        output[10] += 1
        self.assertEqual(y_obs, output)


    def test_observed_distribution_multiple_structures(self):
        train_structures = [
            [
                ('ASP', 'CB', 1, (0, 0, 0)),
                ('LYS', 'CB', 10, (0, 0, 1)),
                ('ARG', 'CB', 10, (0, 10.5, 0))
            ],
            [
                ('ASP', 'CB', 1, (0, 0, 0)),
                ('LYS', 'CB', 10, (0, 0, 1)),
                ('ARG', 'CB', 5, (0, 1.5, 0))
            ],
        ]

        x, y_obs = observed_distribution('ASP', 'ARG', train_structures)
        self.assertEqual(x, list(range(CUTOFF)))

        output = [EPS]* CUTOFF
        output[10] += 1
        output[1] += 1
        self.assertEqual(y_obs, output)

    def test_expected_distribution(self):
        x = [0, 1, 2, 3]
        y_exp = expected_distribution(x)

        correct = [0 + 0 + 1/3, 1 + 1 + 1/3, 4 + 2 + 1/3, 9+3 + 1/3]

        # Ignore any multiplicative factor.
        factor = y_exp[0] / correct[0]
        correct = [x*factor for x in correct]

        self.assertAlmostEqual(y_exp, correct)

    def test_distributions_to_energy(self):
        y_obs = [100, 1000]
        y_exp = [200, 900]
        energy = distributions_to_energy(y_obs, y_exp)
        self.assertAlmostEqual(energy[0], -log(1/2))
        self.assertAlmostEqual(energy[1], -log(10/9))

    def test_scorer_basic(self):
        structure = [
                ('ASP', 'CB', 1, (0, 0, 0)),
                ('ARG', 'CB', 15, (0, 1.5, 0)),
                ('LYS', 'CB', 10, (0, 0, 0.1)),
        ]
        energies = {
            ('ASP', 'ARG'): [10, -1, 0],
            ('ARG', 'ASP'): [10, -1, 0],
            ('ASP', 'LYS'): [8, -2, 0],
            ('LYS', 'ASP'): [8, -2, 0],
            ('ARG', 'LYS'): [9, 0, 0],
            ('LYS', 'ARG'): [9, 0, 0],
        }

        score = scorer(structure, energies)
        self.assertEqual(score, -1 + 8 + 0)

    def test_scorer_over_cutoff(self):
        structure = [
                ('ASP', 'CB', 1, (0, 0, 0)),
                ('ARG', 'CB', 15, (0, 1.5, 0)),
                ('LYS', 'CB', 10, (0, 0, 0.1)),
                ('LYS', 'CB', 20, (0, 0, 100)),
        ]
        energies = {
            ('ASP', 'ARG'): [10, -1, 0],
            ('ARG', 'ASP'): [10, -1, 0],
            ('ASP', 'LYS'): [8, -2, 0],
            ('LYS', 'ASP'): [8, -2, 0],
            ('ARG', 'LYS'): [9, 0, 0],
            ('LYS', 'ARG'): [9, 0, 0],
        }

        score = scorer(structure, energies)

        self.assertEqual(score, -1 + 8 + 0 + 0)

if __name__ == '__main__':
    unittest.main()
