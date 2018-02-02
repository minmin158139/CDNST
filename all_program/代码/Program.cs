using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;

namespace all_program
{
    public class SequenceEle
    {
        public string uid;
        public string tid;
        public int tIndex;
        public DateTime time;
        public int timeIndex;
    }

    class Program
    {
        static void Main(string[] args)
        {
            var inputFile = "../data/insample_music.rpt";
                //读入数据，形成sequence类型的seq

            var seq = File.ReadAllLines(inputFile)
                .Select(dp =>
                {
                    var sp = dp.Split('\t');
                        return new SequenceEle()
                        {
                            uid = sp[0],
                            time = DateTime.Parse(sp[2]),
                            tid = sp[1],


                        };
                   
                }).Where(dp => dp != null)
               .OrderBy(dp => dp.time)
               .ToList();
                
            var inputFile1 = "../data/insample_movie.rpt";

            var seq1 = File.ReadAllLines(inputFile1)
                .Select(dp =>
                {
                    
                
                    var sp = dp.Split('\t');
                        return new SequenceEle()
                        {
                            uid = sp[0],
                            time = DateTime.Parse(sp[2]),
                            tid = sp[1],


                        };
                       }
                       catch (Exception e)
                       {
                           return null;
                       

                }).Where(dp => dp != null)
               .OrderBy(dp => dp.uid).ThenBy(dp => dp.time)
               .ToList();
               
                
                int index = 0;
                var tidIndexDic = seq.Select(dp => dp.tid).Distinct().ToDictionary(dp => dp, dp => index++);

                for (int i = 0; i < seq.Count; i++)
                {
                    seq[i].timeIndex = i;
                    seq[i].tIndex = tidIndexDic[seq[i].tid];



                }

                var selfMatrix = DCN.ComputeSelf(seq);

                int n = seq.Count();
                int m = tidIndexDic.Count();
                int K = 9;



               
             
                int index1 = 0;
                var tidIndexDic1 = seq1.Select(dp => dp.tid).Distinct().ToDictionary(dp => dp, dp => index1++);

                for (int i = 0; i < seq1.Count; i++)
                {
                    seq1[i].timeIndex = i;
                    seq1[i].tIndex = tidIndexDic1[seq1[i].tid];



                }
                var selfMatrix1 = DCN.ComputeSelf(seq1);

                int n1 = seq1.Count();
                int m1 = tidIndexDic1.Count();
                int k1 = 9;
               
                //此时是对此进行训练
              
                for (int kkk = 1; kkk <= 6; kkk++)
                    {
                        Console.WriteLine(kkk + "\t" + "kkkkk");

                        Inference.InitialParameter(m1, n1 - 1, k1);
                        Inference.GibbsSampling_0(seq1.Take(n1 - 1).ToList(), selfMatrix1, 10);

                        Inference.InitialParameter1(m, n, K);
                        Inference.GibbsSampling_0(seq.Take(n - 1).ToList(), selfMatrix, 10);
                       

                        var predicts = new List<Tuple<int, double>>();
                        for (int x = 0; x < m; x++)
                        {

                            double p = 0.0;
                            for (int k = 0; k < K; k++)
                            {
                                try
                                {
                                    p += Math.Pow(1000, 10) * Inference.theta.P[k] * Inference.phi[seq[n - 2].tIndex].P[x] *
                                        Math.Exp(
                                        (-1) * Math.Pow((k + 1 - selfMatrix[n - 2][x]
                                        // / ((double)matrix[i].Max()) * k
                                        ), 2)
                                        );
                                }
                                catch (Exception e)
                                {

                                }
                            }
                            predicts.Add(Tuple.Create(x, p));
                        }
                        predicts = predicts.OrderByDescending(dp => dp.Item2).ToList();


                        double mean = 0.0;
                        for (int k = 0; k < 9; k++)
                        {
                            mean += (k + 1) * Inference.theta.P[k];
                        }
                        int ns = (int)Math.Round(mean);
                        res.Add(seq[0].uid.ToString());
                        res.Add(ns.ToString());




                        var MAP = 0.0;
                        for (int i = 0; i < predicts.Count(); i++)
                        {
                            
                            if (tidIndexDic.Where(dp => dp.Value == predicts.ElementAt(i).Item1).First().Key == seq[n - 1].tid && n > 3)
                            {


                                MAP = i + 1;

                            }

                        }

                        res.Add(MAP.ToString());
                    }
                }
          

        public class DCN
        {
            public static int[][] ComputeSelf(List<SequenceEle> seq)
            {
                int n = seq.Count();
                int m = seq.Select(dp => dp.tid).Distinct().Count();
                int[][] matrix = new int[n][];
                for (int i = 0; i < n; i++)
                    matrix[i] = new int[m];

                for (int i = 0; i < n; i++)
                {
                    if (i == 0)
                    {
                        for (int j = 0; j < m; j++)
                            matrix[i][j] = 1;
                    }
                    else
                    {
                        // 
                        var dic = new Dictionary<int, int>();
                        for (int j = 0; j < m; j++)
                        {
                            dic.Add(j, 1);

                            for (int k = 0; k < i; k++)
                            {
                                var str = seq[k].tid.Split(';');
                                var str1 = seq[j].tid.Split(';');
                                foreach (string str2 in str1)
                                {
                                    if (str.Contains(str2))
                                    {
                                        dic[j]++;
                                    }

                                }
                            }
                            for (int k = 1; k < i; k++)
                            {
                                
                                var str = seq[k].tid.Split(';');
                                var str1 = seq[j].tid.Split(';');
                                var str3 = seq[k - 1].tid.Split(';');
                                var str4 = seq[i - 1].tid.Split(';');
                                bool flag1 = false;
                                bool flag2 = false;
                                foreach (string str2 in str1)
                                {
                                    if (str.Contains(str2))
                                    {
                                        flag1 = true;
                                    }
                                }
                                foreach (string str5 in str4)
                                {
                                    if (str3.Contains(str5))
                                    {
                                        flag2 = true;
                                    }
                                }
                                if (flag1 == true && flag2 == true)
                                    dic[j]++;
                            }
                        }
                        dic = dic.OrderByDescending(dp => dp.Value).ToDictionary(dp => dp.Key, dp => dp.Value);
                        var val = 1;
                        matrix[i][dic.First().Key] = 1;
                        var preFre = dic.First().Value;
                        for (int index = 1; index < dic.Count; index++)
                        {
                            if (dic.ElementAt(index).Value != preFre)
                            {
                                val++;
                                preFre = dic.ElementAt(index).Value;
                            }
                            matrix[i][dic.ElementAt(index).Key] = val;
                        }
                    }
                }

                return matrix;


            }
        }

    }
    public class Inference
    {
        public static int n;
        public static int m;
        public static int k;

        public static double[] alpha;
        public static double[] beta;

        public static MathNet.Numerics.Distributions.Categorical theta;
        public static MathNet.Numerics.Distributions.Categorical[] phi;


        public static int[] zArr;

        public static void InitialParameter(int _m, int _n, int _k)
        {
            m = _m + 1;
            n = _n;
            k = _k;

            alpha = Enumerable.Range(0, m - 1).Select(dp => 1.0).ToArray();
            beta = Enumerable.Range(0, k).Select(dp => 1.0).ToArray();
            MathNet.Numerics.Distributions.Dirichlet dPhi = new MathNet.Numerics.Distributions.Dirichlet(alpha);
            MathNet.Numerics.Distributions.Dirichlet dTheta = new MathNet.Numerics.Distributions.Dirichlet(beta);
            phi = new MathNet.Numerics.Distributions.Categorical[m];
            zArr = new int[n];


            var p = new double[k];
            for (int i = 0; i < k; i++)
                p[i] = 1.0 / k;

            theta = new MathNet.Numerics.Distributions.Categorical(

                p
                );

            for (int i = 0; i < n; i++)
            {
                zArr[i] = theta.Sample();
            }
            for (int i = 0; i < m; i++)
            {
                phi[i] = new MathNet.Numerics.Distributions.Categorical(
                    dPhi.Sample()
                    );
            }
        }



        public static void InitialParameter1(int _m, int _n, int _k)
        {
            m = _m + 1;
            n = _n;
            k = _k;
            alpha = Enumerable.Range(0, m - 1).Select(dp => 1.0).ToArray();
            MathNet.Numerics.Distributions.Dirichlet dPhi = new MathNet.Numerics.Distributions.Dirichlet(alpha);
            phi = new MathNet.Numerics.Distributions.Categorical[m];
            zArr = new int[n];
            for (int i = 0; i < n; i++)
            {
                zArr[i] = theta.Sample();
            }
            for (int i = 0; i < m; i++)
            {
                phi[i] = new MathNet.Numerics.Distributions.Categorical(
                    dPhi.Sample()
                    );
            }


        }

        public static void GibbsSampling_0(List<SequenceEle> seq, int[][] matrix, int iter)
        {
            for (int i = 0; i < iter; i++)
            {
                Console.WriteLine(iter);
                GibbsSamplingEach(seq, matrix);
            }

        }

        public static void GibbsSamplingEach(List<SequenceEle> seq, int[][] matrix)
        {
            for (int i = 0; i < n; i++)
            {
                //draw z
                double[] p = new double[k];
                for (int z = 0; z < k; z++)
                {
                    var theta_tmp = theta.P[z];
                    var phi_tmp = 0.0;
                    if (i == 0)
                        try
                        {
                            phi_tmp = phi[m - 1].P[seq[i].tIndex];
                        }
                        catch (Exception e)
                        {

                        }
                    else
                        try
                        {
                            phi_tmp = phi[seq[i - 1].tIndex].P[seq[i].tIndex];
                        }
                        catch (Exception e)
                        {

                        }
                    if (phi_tmp <= 0)
                        phi_tmp = 1.0 / Math.Pow(10, 100);
                    try
                    {
                        var f = Math.Exp(
                            (-1) * Math.Pow((z + 1 - matrix[i][seq[i].tIndex]
                            ), 2)
                            );
                        p[z] = Math.Pow(10, 10) * theta_tmp * phi_tmp * f;
                    }
                    catch (Exception e)
                    {

                    }
                }
                try
                {
                    MathNet.Numerics.Distributions.Categorical dis = new MathNet.Numerics.Distributions.Categorical(p);
                    zArr[i] = dis.Sample();
                }
                catch (Exception e)
                {

                }

                var beta_tmp = new double[k];
                beta.CopyTo(beta_tmp, 0);
                foreach (var g in zArr.GroupBy(dp => dp))
                {
                    beta_tmp[g.Key] += g.Count();
                }
                var d_tmp = new MathNet.Numerics.Distributions.Dirichlet(beta_tmp);
                theta = new MathNet.Numerics.Distributions.Categorical(d_tmp.Sample());

                for (int j = 0; j < m; j++)
                {
                    p = Enumerable.Range(0, m).Select(dp => 1.0).ToArray();
                    for (int ii = 0; ii < n; ii++)
                    {

                        if (j == m - 1 && ii == 0)
                        {
                            try
                            {


                                var z = zArr[ii];

                                var f = Math.Exp(
                        (-1) * Math.Pow((z - matrix[ii][seq[ii].tIndex] / ((double)matrix[ii].Max()) * k), 2)
                        );

                                p[seq[ii].tIndex] *= phi[j].P[seq[ii].tIndex] * f;

                            }
                            catch (Exception e)
                            {

                            }
                        }


                        else
                        {

                            if (ii == 0)
                                continue;
                            try
                            {
                                if (seq[ii - 1].tIndex == j)
                                {
                                    var z = zArr[ii];

                                    var f = Math.Exp(
                            (-1) * Math.Pow((z + 1 - matrix[ii][seq[ii].tIndex]
                            ), 2)
                            );

                                    p[seq[ii].tIndex] *= phi[j].P[seq[ii].tIndex] * f;
                                }
                            }
                            catch (Exception e)
                            {

                            }
                        }
                    }

                    for (int s = 0; s < p.Length; s++)
                        if (p[s] == 1.0 || p[s] <= 0)
                            p[s] = 1.0 / Math.Pow(10, 100);
                    phi[j] = new MathNet.Numerics.Distributions.Categorical(p);
                }
            }
        }





