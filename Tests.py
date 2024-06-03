from unittest import main as startTests, TestCase
from Vec import *
from Poly import *
from Mat import *

class Test(TestCase):
    def testPoly(self):
        P = Poly(coefs=[1,2,3]) 
        Q = Poly(coefs=[1,-1,2,-2,0])
        self.assertEqual(P.deg(), 2)
        self.assertEqual(Q.deg(), 3)
        self.assertEqual(P + Q, Poly(coefs=[2,1,5,-2]))
        self.assertEqual(P - Q, Poly(coefs=[0,3,1,2]))
        self.assertEqual(P * Q, Poly(coefs=[1,1,3,-1,2,-6]))
        self.assertEqual(Q * P, Poly(coefs=[1,1,3,-1,2,-6]))
        (D,R) = P.euclidian_div(Q)
        self.assertAlmostEqual(abs(P * D + R - Q), 0)
        self.assertEqual(P(1), 6)
        self.assertEqual(P(2), 17)
        self.assertEqual(P(Idt(2)), 6 * Idt(2))
        self.assertEqual(P(2 * Idt(2)), 17 * Idt(2))

    def testMat1(self):
        M = Idt(3)
        for i in range(3):
            for j in range(3):
                x = 1 if i == j else 0
                self.assertEqual(M[i][j], x)
                self.assertEqual(M.lines[i][j], x)
                self.assertEqual(M.cols[j][i], x)

    def testMat2(self):
        M = Idt(3)
        self.assertTrue(M == (M + M) - M)
        self.assertTrue(M == M * 2 - M)
        self.assertTrue(M == M*M)

    def testMatInv(self):
        M = Mat(coefs=[[2,1],[-1,0]])
        self.assertAlmostEqual(abs(M * M.inv() - Idt(2)), 0)
        self.assertAlmostEqual(abs(M.inv() * M - Idt(2)), 0)
        self.assertEqual(M.carPoly()(M), Zero(2))
        M = Mat(coefs=[[1,2,3],[4,5,6],[7,8,1]])
        self.assertAlmostEqual(abs(M * M.inv() - Idt(3)), 0)
        self.assertAlmostEqual(abs(M.inv() * M - Idt(3)), 0)
        self.assertEqual(M.carPoly()(M), Zero(3))

startTests()
    
