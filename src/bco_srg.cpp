#include "../include/bco_srg.hpp"
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sys/time.h>

using namespace std;

/*  normalization function to convert randomly
    generated number into an integer within
    given range*/
extern int normal(int, int, int);

/*
    function for getting random element from an
    unordered_set in constant time
*/
extern int random_element_set(unordered_set<int> &, int);

/*  function for calculating spectrum of a graph
    using jacobi algorithm
*/
extern void jacobi_eigenvalue(int, double[], int, double[], double[], int *, int *);

/* function for calculating euclidean distance */
extern double spectral_distance(double[], double[], int);

/* function for printing out the solution */
void solution(double *As, double *sp, int n, double time_taken, double ygmin, int eval) {
    fprintf(matr_spectrum, "Matrix of adjacency: \n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(matr_spectrum, "%lf ", As[i * n + j]);
        }
        fprintf(matr_spectrum, "\n");
    }
    fprintf(matr_spectrum, "Spectrum: \n");
    for (int i = 0; i < n; i++) {
        fprintf(matr_spectrum, "%lf ", sp[i]);
    }
    fprintf(matr_spectrum, "\n");

    fprintf(fout, "%.4lf\t%.4lf\t%d", ygmin, time_taken, eval);
}

/* MAIN */
int bco_srg(int seed, int bees, int NC) {
    if ((fin = fopen("res/input.txt", "r")) == NULL) {
        printf("Open INPUT.DAT.\n");
        exit(0);
    }

    if ((fout = fopen("res/output.txt", "a")) == NULL) {
        printf("Open %s.\n", "output.dat");
        exit(0);
    }
    if ((matr_spectrum = fopen("res/graphs_output.txt", "a")) == NULL) {
        printf("Open %s.\n", "graphs_output.dat");
        exit(0);
    }
    fprintf(fout, "\n%d\t%d\t%d\t%d\t%s\t", bees, NC, 537, seed, "10_2");
    int best_bee;   /* index of bee that generated best (partial) solution */
    int R;          /* number of recruiters (loyal bees) */
    int pass_count; /* forward pass counter */

    int loy[MAX_NO_BEES]; /* loyalty indicator array */
    int non_loy[MAX_NO_BEES];

    int aux, aux1, aux2, loc; /* auxiliary variables */
    long o, o_range;          /* number of edges that will be replaced */

    int stop_ind; /* stopping criteria indicator: 0-no. iter; 1-no. unimp. iter; 2-runtime */
    long num_iter;
    double run_time; /* allowed CPU runtime in seconds */

    long n_it,       /* current iteration number */
        improvement, /* indicator of the current best solution improvement */
        n_nonimp_it; /* current number of iterations without improvement */
    long terminate;  /* stopping criteria indicator */

    double max_y, min_y; /* max and min partial solution (among all bees in one forward pass) */

    long value, bingo;     /* integer roulet variables */
    double dvalue, dbingo; /* probability roulet variables */
    double maxV_norm;      /* maximum of normalized values for bees partial solutions */
    double sum_V_loy;      /* sum of normalized values for bees partial solutions */

    double V_norm[MAX_NO_BEES];   /* normalized values for bees partial solutions */
    double p_loy[MAX_NO_BEES];    /* probability that bee is loyal */
    double Droulete[MAX_NO_BEES]; /* probability that recruiter is chosen */

    struct timeval start, end;

    int it_max = 100000; /* max no of iterations in jacobi algorithm */
	int eval = 0;	/* number of evaluations of objective function*/
    long fv, sv, sv_ind; /* local variables for selecting first and second vertex when generating edge*/
    long fv_out, sv_out, fv_in, sv_in;

    // Time variables.
    clock_t t, t1, t2;
    double min_t;
    stop_ind = 0;
    num_iter = 537;
    fscanf(fin, "%d", &n);
    for (int i = 0; i < n; i++) {
        fscanf(fin, "%lf", &C[i]);
    }

    double m_pom = 0;
    for (int i = 0; i < n; i++) {
        m_pom += C[i] * C[i];
    }
    m = (int)round((m_pom / 2));
    n_it = n_nonimp_it = 0;
    terminate = 0;

    gettimeofday(&start, NULL);
    t1 = clock();
    bool first_it = true;
    bool first_init = true;
    ygmin = DBL_MAX;

    // variables to be passed as parameters for jacobi_eigenvalue function
    double *v = (double *)malloc(sizeof(double) * (n * n));
    if (v == NULL) {
        printf("Malloc");
        return 1;
    }
    double *sp = (double *)malloc(sizeof(double) * n);
    if (sp == NULL) {
        printf("malloc");
        return 1;
    }
    int it_num;
    int rot_num;

    srand(seed);
    while (!terminate) { /* until stopping criteria */
        max_y = (double)0.0;
        min_y = DBL_MAX;
        n_it++; /* new iteration starts */

        // neighbors and non neighbors sets for each bee
        vector<vector<unordered_set<int>>> neigh_sets(MAX_NO_BEES, vector<unordered_set<int>>(MAX_NO_NODES));
        vector<vector<unordered_set<int>>> non_neigh_sets(MAX_NO_BEES, vector<unordered_set<int>>(MAX_NO_NODES));

        // sets of non isolated nodes
        vector<unordered_set<int>> non_iso(MAX_NO_BEES);

        // set of nodes with less than n-1 edges
        vector<unordered_set<int>> not_compl_nodes(MAX_NO_BEES);
        /* construction of initial solutions */
        for (int b = 0; b < bees; b++) {
            for (int i = 0; i < n; i++) {
                not_compl_nodes[b].insert(i);
                for (int j = i; j < n; j++) {
                    As[b][i * n + j] = 0.0;
                    As[b][j * n + i] = 0.0;
                    if (i != j) {
                        non_neigh_sets[b][i].insert(j);
                        non_neigh_sets[b][j].insert(i);
                    }
                }
            }

            for (int e = 0; e < m; e++) {

                // selecting random edge
                aux1 = rand();
                aux1 = normal(aux1, 0, n - 1); // first vertex
                fv = random_element_set(not_compl_nodes[b], aux1);
                aux2 = rand();
                aux2 = normal(aux2, 0, n - 1);
                sv = random_element_set(non_neigh_sets[b][fv], aux2);
                // updating data
                As[b][fv * n + sv] = 1.0;
                As[b][sv * n + fv] = 1.0;
                neigh_sets[b][fv].insert(sv);
                non_neigh_sets[b][fv].erase(sv);
                neigh_sets[b][sv].insert(fv);
                non_neigh_sets[b][sv].erase(fv);
                if (neigh_sets[b][fv].size() == (n - 1)) {
                    not_compl_nodes[b].erase(fv);
                }
                if (neigh_sets[b][sv].size() == (n - 1)) {
                    not_compl_nodes[b].erase(sv);
                }
            }

            for (int i = 0; i < n; i++) {
                if (neigh_sets[b][i].size() > 0) {
                    non_iso[b].insert(i);
                }
            }

            jacobi_eigenvalue(n, As[b], it_max, v, sp, &it_num, &rot_num);

            y[b] = spectral_distance(C, sp, n);
            if (b == 0) {

                max_y = y[0];
                min_y = y[0];
                best_bee = 0;
                if (first_it || (ygmin - y[0]) > EPS) {
                    ygmin = y[0];
                    t2 = clock();
                    min_t = (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
                    if (ygmin < EPS) {
                        solution(As[b], sp, n, min_t, ygmin, eval);
                        free(v);
                        free(sp);
                        if (fclose(fout) == EOF) {
                            printf("fclose");
                        }
                        if (fclose(fin) == EOF) {
                            printf("fclose");
                        }
                        return 0;
                    }
                    improvement = 1;
                    for (int i = 0; i < n; i++) {
                        for (int j = i; j < n; j++) {
                            Agmin[i * n + j] = As[best_bee][i * n + j];
                            Agmin[j * n + i] = As[best_bee][j * n + i];
                        }
                    }

                    first_it = false;
                }
            } else {
                if (max_y < y[b]) max_y = y[b];
                if (min_y > y[b]) {
                    best_bee = b;
                    min_y = y[b];
                }
                if ((ygmin - y[best_bee]) > EPS) {
                    t2 = clock();
                    min_t = (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
                    ygmin = y[best_bee];
                    if (ygmin < EPS) {
                        solution(As[b], sp, n, min_t, ygmin, eval);
                        free(v);
                        free(sp);
                        fclose(fout);
                        fclose(fin);
                        return 0;
                    }
                    improvement = 1;
                    for (int i = 0; i < n; i++) {
                        for (int j = i; j < n; j++) {
                            Agmin[i * n + j] = As[best_bee][i * n + j];
                            Agmin[j * n + i] = As[best_bee][j * n + i];
                        }
                    }
                }
            }

        } /* end of construction phase */
		eval += bees;
        for (pass_count = 1; pass_count <= NC; pass_count++) {
            /* BACKWARD PASS */
            if (max_y == min_y) {
                for (int b = 0; b < bees; b++) {
                    loy[b] = b;
                }
            } else {
                dvalue = (double)(max_y - min_y);
                for (int b = 0; b < bees; b++) {
                    V_norm[b] = (double)(max_y - y[b]) / dvalue;
                }
                // counters for loyal and non-loyal bees
                int no_loy = 0;
                int no_non_loy = 0;
                sum_V_loy = 0;
                for (int b = 0; b < bees; b++) {
                    p_loy[b] = V_norm[b];
                    aux = rand();
                    dbingo = (double)aux / (double)RAND_MAX;
                    if (dbingo < p_loy[b]) {
                        loy[no_loy++] = b;
                        sum_V_loy += V_norm[b];
                    } else {
                        non_loy[no_non_loy++] = b;
                    }
                }

                for (int b = 0; b < no_non_loy; b++) {
                    aux = rand(); /* for each non-loyal bee */
                    dbingo = (double)aux / (double)RAND_MAX;

                    // selecting recruiter
                    int bb = 0;

                    double dvalue = V_norm[loy[0]] / sum_V_loy;
                    while (bb < no_loy && dbingo - dvalue > EPS) {
                        dvalue += V_norm[loy[++bb]] / sum_V_loy;
                    }
                    int recruiter = loy[bb];
                    int follower = non_loy[b];
                    /* copying solution from recruiter to follower */
                    y[follower] = y[recruiter];
                    for (int i = 0; i < n; i++) {

                        neigh_sets[follower][i] = neigh_sets[recruiter][i];
                        non_neigh_sets[follower][i] = non_neigh_sets[recruiter][i];
                        non_iso[follower] = non_iso[recruiter];
                        not_compl_nodes[follower] = not_compl_nodes[recruiter];

                        for (int j = i; j < n; j++) {
                            As[follower][(i * n) + j] = As[recruiter][(i * n) + j];
                            As[follower][(j * n) + i] = As[recruiter][(j * n) + i];
                        }
                    }
                }
            }
            /* END OF BACKWARD PASS */

            for (int b = 0; b < bees; b++) { /* FORWARD PASS */
                int o_range = m / 4;
                aux = rand();
                int o = normal(aux, 1, o_range);
                for (int i = 0; i < o; i++) {
                    aux1 = rand();
                    aux1 = normal(aux1, 0, n - 1);
                    fv_out = random_element_set(non_iso[b], aux1); // first vertex of edge to be erased
                    aux2 = rand();
                    aux2 = normal(aux2, 0, n - 1);
                    sv_out = random_element_set(neigh_sets[b][fv_out], aux2); // second vertex of edge to be erased
                    aux1 = rand();
                    aux1 = normal(aux1, 0, n - 1);
                    fv_in = random_element_set(not_compl_nodes[b], aux1); // first vertex of edge to be added
                    aux2 = rand();
                    aux2 = normal(aux2, 0, n - 1);
                    sv_in = random_element_set(non_neigh_sets[b][fv_in], aux2); // second vertex of edge to be added

                    // adding edge
                    As[b][fv_in * n + sv_in] = 1.0;
                    As[b][sv_in * n + fv_in] = 1.0;
                    neigh_sets[b][fv_in].insert(sv_in);
                    non_neigh_sets[b][fv_in].erase(sv_in);
                    neigh_sets[b][sv_in].insert(fv_in);
                    non_neigh_sets[b][sv_in].erase(fv_in);
                    if (neigh_sets[b][fv_in].size() == 1) {
                        non_iso[b].insert(fv_in);
                    }
                    if (neigh_sets[b][sv_in].size() == 1) {
                        non_iso[b].insert(sv_in);
                    }
                    if (neigh_sets[b][fv_in].size() == (n - 1)) {
                        not_compl_nodes[b].erase(fv_in);
                    }
                    if (neigh_sets[b][sv_in].size() == (n - 1)) {
                        not_compl_nodes[b].erase(sv_in);
                    }

                    // removing edge
                    As[b][fv_out * n + sv_out] = 0.0;
                    As[b][sv_out * n + fv_out] = 0.0;
                    neigh_sets[b][fv_out].erase(sv_out);
                    non_neigh_sets[b][fv_out].insert(sv_out);
                    neigh_sets[b][sv_out].erase(fv_out);
                    non_neigh_sets[b][sv_out].insert(fv_out);
                    if (neigh_sets[b][fv_out].size() == (n - 2)) {
                        not_compl_nodes[b].insert(fv_out);
                    }
                    if (neigh_sets[b][sv_out].size() == (n - 2)) {
                        not_compl_nodes[b].insert(sv_out);
                    }
                    if (neigh_sets[b][fv_out].size() == 0) {
                        non_iso[b].erase(fv_out);
                    }
                    if (neigh_sets[b][sv_out].size() == 0) {
                        non_iso[b].erase(sv_out);
                    }
                }
                // calculating objective function
                jacobi_eigenvalue(n, As[b], it_max, v, sp, &it_num, &rot_num);
                y[b] = spectral_distance(C, sp, n);

                if (b == 0) {
                    max_y = y[0];
                    min_y = y[0];
                    best_bee = 0;
                }
                if (max_y < y[b]) max_y = y[b];
                if (min_y > y[b]) {
                    best_bee = b;
                    min_y = y[b];
                }
                if ((ygmin - y[best_bee]) > EPS) {
                    ygmin = y[best_bee];
                    improvement = 1;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            Agmin[i * n + j] = As[best_bee][i * n + j];
                        }
                    }
                    t2 = clock();
                    min_t = (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
                    if (ygmin < EPS) {
                    	eval = eval + (pass_count - 1) * bees + b;
                        solution(As[b], sp, n, min_t, ygmin, eval);
                        free(v);
                        free(sp);
                        if (fclose(fout) == EOF) {
                            printf("fclose");
                        }
                        if (fclose(fin) == EOF) {
                            printf("fclose");
                        }

                        return 0;
                    }
                }

            } /* END OF FORWARD PASS */
        }     /* NC */
        eval += bees*NC;	

        if (improvement)
            n_nonimp_it = 0;
        else {
            n_nonimp_it++;
        }
        /* checking if stopping criteria is satisfied */
        switch (stop_ind) {
        case 0:
            if (n_it >= num_iter) {
                terminate = 1;
                t2 = clock();
                min_t = (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
                jacobi_eigenvalue(n, Agmin, it_max, v, sp, &it_num, &rot_num);
                solution(Agmin, sp, n, min_t, spectral_distance(C, sp, n), eval);
                fflush(fout);
                fflush(matr_spectrum);
                fclose(fout);
                fclose(fin);
                free(v);
                free(sp);
            }
            break;
        case 1:
            if (n_nonimp_it >= num_iter) {
                fprintf(fout, "Failed to find solution (max number of iterations without improvement)");
                terminate = 1;
            }
            break;
        case 2:
            t2 = clock();
            t = (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
            if (t > run_time) {
                terminate = 1;
                fprintf(fout, "Failed to find solution (time) t=%lf", t);
            }
            break;
        }
    } /* first while loop */

    return 0;
}
