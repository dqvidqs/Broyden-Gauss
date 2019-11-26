using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
//----------14 VARIANTAS-----------
namespace lab_02
{
    public partial class Form1 : Form
    {
        List<Timer> Timerlist = new List<Timer>();

        public Form1()
        {
            InitializeComponent();
            Initialize();
        }
        // ---------------------------------------------- GAUSO METODAS ----------------------------------------------
        double[,] A = { { 2, 5, 1, 2 }, { -2, 0, 3, 5 }, { 1, 0, -1, 1 }, { 0, 5, 4, 7 } }; // koeficientu matrica
        double[] B = { 14, 10, 4, 24 };// koeficientu vektorius
        /*double[,] A = { { 1, 1, 1, 1 }, { 1, -1, -1, 1 }, { 2, 1, -1, 2 }, { 3, 1, 2, -1 } }; // koeficientu matrica
        double[] B = { 2, 0, 9, 7 };// koeficientu vektorius*/
        private void button2_Click(object sender, EventArgs e)
        {
            ClearForm();
            Gaussmethod();
        }

        private void Gaussmethod()
        {
            double[,] AB = new double[B.Length, A.GetLength(0) + 1];
            //sukuria nauja matrica
            for (int i = 0; i < AB.GetLength(0); i++)
            {
                for (int j = 0; j < AB.GetLength(1); j++)
                {

                    if (j >= A.GetLength(0))
                    {
                        AB[i, j] = B[i];
                    }
                    else
                    {
                        AB[i, j] = A[i, j];
                    }
                }
            }
            PrintToRichBoxGauseMethod(AB, -1, -1);
            //Gauso metodas
            for (int i = 0; i < AB.GetLength(0) - 1; i++)
            {
                for (int j = i + 1; j < AB.GetLength(0); j++)
                {
                    double X = AB[j, i] / AB[i, i];
                    int Index = j;
                    for (int k = i; k < AB.GetLength(1); k++)
                    {
                        AB[j, k] = AB[j, k] - AB[i, k] * X;
                    }
                    PrintToRichBoxGauseMethod(AB, X, Index);
                }
            }
            //x gaunami
            double[] x = new double[AB.GetLength(0)];
            for (int i = AB.GetLength(0) - 1; i > -1; i--)
            {
                double sum = 0;
                for (int j = AB.GetLength(1) - 2; j > i; j--)
                {
                    sum = sum - (x[j] * AB[i, j]);
                }
                if (AB[i, i] == 0 && AB[i, AB.GetLength(1) - 1] == 0)//laisvai parenkamas narys
                {
                    x[i] = 0;
                }
                else
                {
                    x[i] = (AB[i, AB.GetLength(1) - 1] + sum) / AB[i, i];
                }
            }
            //x atspausdinami
            for (int i = 0; i < x.Length; i++)
            {
                richTextBox1.AppendText(string.Format("x{1} = {0}\n", x[i], i+1));
            }
            richTextBox1.AppendText("\n");
            //patikrinimas
            richTextBox1.AppendText("Patikrinimas (liekana)\n");
            for (int i = 0; i < A.GetLength(1); i++)
            {
                double sum = 0;
                for (int j = 0; j < x.Length; j++)
                {
                    sum = sum + x[j] * A[i, j];
                }
                sum = sum - B[i];
                richTextBox1.AppendText(string.Format("{0} eilute ats: {1}\n",i+1,sum));
            }
            richTextBox1.AppendText("\n");
        }
        private void PrintToRichBoxGauseMethod(double[,] AB, double X, int Index)
        {
            richTextBox1.AppendText(string.Format("Matrix {0}x{1}\n", AB.GetLength(0), AB.GetLength(1)));
            for (int i = 0; i < AB.GetLength(0); i++)
            {
                for (int j = 0; j < AB.GetLength(1); j++)
                {
                    richTextBox1.AppendText(string.Format("{0,-4} ",AB[i, j]));
                    if (Index == i && j + 1 == AB.GetLength(1))
                    {
                        richTextBox1.AppendText(string.Format("| /{0}",X));
                    }
                }
                richTextBox1.AppendText("\n");
            }
            richTextBox1.AppendText("\n");
        }
        // ----------------------------------------------  BROIDENO METODAS ----------------------------------------------
        private int Iter_max = 100;
        private int Iter = 0;
        private double[] x = { -3.5, -8 };
        private double[,] A_br = { { 0, 0 }, { 0, 0 } };
        private double[] ff;
        Matrix<double> B_ARTINIAI;
        Matrix<double> FF_TEMP;
        double eps = 1e-6;

        private double[] Funk(double[] x_temp)
        {
            double[] ff = { 0, 0 };
            ff[0] = x_temp[0] * (x_temp[1] + 2 * Math.Cos(x_temp[0])) - 1;
            ff[1] = Math.Pow(x_temp[0], 4) + Math.Pow(x_temp[1], 4)-64;
            return ff;
        }

        private void button3_Click(object sender, EventArgs e)
        {
            BroidenMethod();
        }
        private void BroidenMethod()
        {
            //--
            double DX = (Math.Abs(x[0]) + Math.Abs(x[1])) / 100000;
            double[] f0 = Funk(x);

            for (int i = 0; i < x.Length; i++)
            {
                double[] x_temp = { x[0], x[1] };
                x_temp[i] = x_temp[i] + DX;
                double[] f1 = Funk(x_temp);
                for (int j = 0; j < f1.Length; j++)
                {
                    A_br[i, j] = (f1[j] - f0[j]) / DX;
                }
            }
            //--
            ff = Funk(x);
            double[,] temp_ff = { { ff[0] }, { ff[1] } };
            B_ARTINIAI = Matrix<double>.Build.DenseOfArray(A_br);
            FF_TEMP = Matrix<double>.Build.DenseOfArray(temp_ff);
            timer1.Enabled = true;
            timer1.Interval = 50; // timer1 intervalas milisekundemis
            timer1.Start();
        }
        private void timer1_Tick(object sender, EventArgs e)
        {
            if (Iter < Iter_max)
            {
                Iter++;
                //laikinas ff
                double[,] temp_ff = { { ff[0] }, { ff[1] } };
                FF_TEMP = Matrix<double>.Build.DenseOfArray(temp_ff); ;

                var deltax = B_ARTINIAI.Solve(FF_TEMP);

                double[] x1 = {0,0 };
                for (int i = 0; i < 2; i++)
                {
                    x1[i] = x[i] + deltax.At(i, 0);
                }
                double[] ff1 = Funk(x1);
                double[,] temp_ff1 = { { ff1[0] }, { ff1[1] } };

                Matrix<double> FF1_TEMP = Matrix<double>.Build.DenseOfArray(temp_ff1);

                var temp_up = B_ARTINIAI + (FF1_TEMP - FF_TEMP - B_ARTINIAI * deltax) * deltax.Transpose();
                var temp_down = deltax.Transpose() * deltax;
                double temp_down_0 = temp_down.At(0, 0);
                B_ARTINIAI = temp_up * temp_down_0;

                double[,] temp_x = { { x[0] }, { x[1] } };
                Matrix<double> X_TEMP = Matrix<double>.Build.DenseOfArray(temp_x);
                var tikslumas = deltax.InfinityNorm() / X_TEMP.InfinityNorm() + deltax.InfinityNorm();
                if (tikslumas < eps)
                { timer1.Stop(); }
                richTextBox1.AppendText(string.Format("{0} {1}\n",Iter,tikslumas));
                //--
                ff = ff1;
                x = x1;
            }
            else
            {
                richTextBox1.AppendText("Skaičiavimai baigti\n");
                timer1.Stop();
            }

        }
        // ---------------------------------------------- KITI METODAI ----------------------------------------------
        /// <summary>
        /// Uždaroma programa
        /// </summary>
        private void button1_Click(object sender, EventArgs e)
        {
            Close();
        }
        
        /// <summary>
        /// Išvalomas grafikas ir consolė
        /// </summary>
        private void button4_Click(object sender, EventArgs e)
        {
            ClearForm();
        }
        
        public void ClearForm()
        {
            richTextBox1.Clear(); // isvalomas richTextBox1
            // sustabdomi timeriai jei tokiu yra
            foreach (var timer in Timerlist)
            {
                timer.Stop();
            }

            // isvalomos visos nubreztos kreives
            chart1.Series.Clear();
        }

    }
}
