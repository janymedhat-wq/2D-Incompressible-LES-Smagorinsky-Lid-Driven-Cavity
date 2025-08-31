// 2D Incompressible LES (Smagorinsky) – Lid-Driven Cavity
// les_cavity.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IDX(i,j,Nx) ((i) + (Nx)*(j))

typedef struct {
    int Nx, Ny;
    double Lx, Ly, dx, dy, dt, Re, rho, Cs, Ulid;
    int    it_max, poisson_it, save_every;
    double beta_poisson; // SOR factor ~1.7
} Params;

static void zero(double *a, int n){ for(int i=0;i<n;i++) a[i]=0.0; }
static void copy(double *dst, const double *src, int n){ for(int i=0;i<n;i++) dst[i]=src[i]; }

int main(){
    Params P = {
        .Nx=128, .Ny=128, .Lx=1.0, .Ly=1.0,
        .dt=5e-4, .Re=10000.0, .rho=1.0, .Cs=0.14, .Ulid=1.0,
        .it_max=20000, .poisson_it=200, .save_every=1000, .beta_poisson=1.7
    };
    P.dx = P.Lx/(P.Nx-1); P.dy = P.Ly/(P.Ny-1);
    double nu = P.Ulid * P.Lx / P.Re;

    int N = P.Nx * P.Ny;
    double *u = malloc(sizeof(double)*N);
    double *v = malloc(sizeof(double)*N);
    double *p = malloc(sizeof(double)*N);
    double *u_star = malloc(sizeof(double)*N);
    double *v_star = malloc(sizeof(double)*N);
    double *rhs = malloc(sizeof(double)*N);
    double *nu_t = malloc(sizeof(double)*N);

    zero(u,N); zero(v,N); zero(p,N);
    zero(u_star,N); zero(v_star,N); zero(rhs,N); zero(nu_t,N);

    const double dx = P.dx, dy = P.dy, dt = P.dt;
    const double Cs = P.Cs, Ulid=P.Ulid, rho=P.rho;

    FILE *fp = NULL;

    for(int it=0; it<P.it_max; ++it){
        // 1) Boundary conditions (no-slip; moving lid at j=Ny-1)
        for(int i=0;i<P.Nx;i++){
            u[IDX(i,0,P.Nx)] = 0.0; v[IDX(i,0,P.Nx)] = 0.0;
            u[IDX(i,P.Ny-1,P.Nx)] = Ulid; v[IDX(i,P.Ny-1,P.Nx)] = 0.0;
        }
        for(int j=0;j<P.Ny;j++){
            u[IDX(0,j,P.Nx)] = 0.0; v[IDX(0,j,P.Nx)] = 0.0;
            u[IDX(P.Nx-1,j,P.Nx)] = 0.0; v[IDX(P.Nx-1,j,P.Nx)] = 0.0;
        }

        // 2) Compute Smagorinsky eddy viscosity nu_t
        double Delta = dx; // uniform grid
        double Cs2Delta2 = Cs*Cs*Delta*Delta;
        for(int j=1;j<P.Ny-1;j++){
            for(int i=1;i<P.Nx-1;i++){
                int k = IDX(i,j,P.Nx);
                double du_dx = (u[IDX(i+1,j,P.Nx)] - u[IDX(i-1,j,P.Nx)])/(2*dx);
                double dv_dy = (v[IDX(i,j+1,P.Nx)] - v[IDX(i,j-1,P.Nx)])/(2*dy);
                double du_dy = (u[IDX(i,j+1,P.Nx)] - u[IDX(i,j-1,P.Nx)])/(2*dy);
                double dv_dx = (v[IDX(i+1,j,P.Nx)] - v[IDX(i-1,j,P.Nx)])/(2*dx);

                double Sxx = du_dx;
                double Syy = dv_dy;
                double Sxy = 0.5*(du_dy + dv_dx);
                double Smag = sqrt(2.0*(Sxx*Sxx + Syy*Syy + 2.0*Sxy*Sxy));
                nu_t[k] = Cs2Delta2 * Smag;
            }
        }

        // 3) Build intermediate velocities (explicit advection, viscous with nu+nu_t)
        for(int j=1;j<P.Ny-1;j++){
            fDX(i,j+1,P.Nx)], vS=v[IDX(i,j-1,P.Nx)];

                double du_dx = (uE - uW)/(2*dx);
                double du_dy = (uN - uS)/(2*dy);
                double dv_dx = (vE - vW)/(2*dx);
                double dv_dy = (vN - vS)/(2*dy);

                double adv_u = u[k]*du_dx + v[k]*du_dy;
                double adv_v = u[k]*dv_dx + v[k]*dv_dy;

                // variable viscosity diffusion: ∇·[(nu+nu_t)∇u]
                double nue = 0.5*(nu+nu_t[k] + nu+nu_t[IDX(i+1,j,P.Nx)]);
                double nuw = 0.5*(nu+nu_t[k] + nu+nu_t[IDX(i-1,j,P.Nx)]);
                double nun = 0.5*(nu+nu_t[k] + nu+nu_t[IDX(i,j+1,P.Nx)]);
                double nus = 0.5*(nu+nu_t[k] + nu+nu_t[IDX(i,j-1,P.Nx)]);

                double lap_u = (nue*(uE - u[k]) - nuw*(u[k]-uW))/(dx*dx)
                             + (nun*(uN - u[k]) - nus*(u[k]-uS))/(dy*dy);
                double lap_v = (nue*(vE - v[k]) - nuw*(v[k]-vW))/(dx*dx)
                             + (nun*(vN - v[k]) - nus*(v[k]-vS))/(dy*dy);

                u_star[k] = u[k] + dt * (-adv_u + lap_u);
                v_star[k] = v[k] + dt * (-adv_v + lap_v);
            }
        }

        // Re-apply BCs to star fields
        for(int i=0;i<P.Nx;i++){
            u_star[IDX(i,0,P.Nx)] = 0.0; v_star[IDX(i,0,P.Nx)] = 0.0;
            u_star[IDX(i,P.Ny-1,P.Nx)] = Ulid; v_star[IDX(i,P.Ny-1,P.Nx)] = 0.0;
        }
        for(int j=0;j<P.Ny;j++){
            u_star[IDX(0,j,P.Nx)] = 0.0; v_star[IDX(0,j,P.Nx)] = 0.0;
            u_star[IDX(P.Nx-1,j,P.Nx)] = 0.0; v_star[IDX(P.Nx-1,j,P.Nx)] = 0.0;
        }

        // 4) Pressure Poisson RHS = rho/dt * div(u*)
        for(int j=1;j<P.Ny-1;j++){
            for(int i=1;i<P.Nx-1;i++){
                int k = IDX(i,j,P.Nx);
                double du_dx = (u_star[IDX(i+1,j,P.Nx)] - u_star[IDX(i-1,j,P.Nx)])/(2*dx);
                double dv_dy = (v_star[IDX(i,j+1,P.Nx)] - v_star[IDX(i,j-1,P.Nx)])/(2*dy);
                rhs[k] = rho/dt * (du_dx + dv_dy);
            }
        }

        // 5) Solve Poisson: ∇^2 p = rhs (SOR)
        for(int itp=0; itp<P.poisson_it; ++itp){
            for(int j=1;j<P.Ny-1;j++){
                for(int i=1;i<P.Nx-1;i++){
                    int k = IDX(i,j,P.Nx);
                    double p_new = ((p[IDX(i+1,j,P.Nx)] + p[IDX(i-1,j,P.Nx)])/(dx*dx) +
                                    (p[IDX(i,j+1,P.Nx)] + p[IDX(i,j-1,P.Nx)])/(dy*dy) - rhs[k]) /
                                   (2.0*(1.0/(dx*dx) + 1.0/(dy*dy)));
                    p[k] = (1.0 - P.beta_poisson)*p[k] + P.beta_poisson*p_new;
                }
            }
            // Neumann BC: dp/dn = 0 -> copy interior
            for(int i=0;i<P.Nx;i++){ p[IDX(i,0,P.Nx)] = p[IDX(i,1,P.Nx)]; p[IDX(i,P.Ny-1,P.Nx)] = p[IDX(i,P.Ny-2,P.Nx)]; }
            for(int j=0;j<P.Ny;j++){ p[IDX(0,j,P.Nx)] = p[IDX(1,j,P.Nx)]; p[IDX(P.Nx-1,j,P.Nx)] = p[IDX(P.Nx-2,j,P.Nx)]; }
        }

        // 6) Projection: u = u* - dt/rho * grad p
        for(int j=1;j<P.Ny-1;j++){
            for(int i=1;i<P.Nx-1;i++){
                int k = IDX(i,j,P.Nx);
                double dp_dx = (p[IDX(i+1,j,P.Nx)] - p[IDX(i-1,j,P.Nx)])/(2*dx);
                double dp_dy = (p[IDX(i,j+1,P.Nx)] - p[IDX(i,j-1,P.Nx)])/(2*dy);
                u[k] = u_star[k] - dt/rho * dp_dx;
                v[k] = v_star[k] - dt/rho * dp_dy;
            }
        }

        // 7) Output occasionally (CSV for quick plots)
        if(it % P.save_every == 0){
            char name[128];
            snprintf(name,sizeof(name),"fields_%06d.csv", it);
            fp = fopen(name,"w");
            fprintf(fp,"# x,y,u,v,p\n");
            for(int j=0;j<P.Ny;j++){
                for(int i=0;i<P.Nx;i++){
                    double x = i*dx, y = j*dy;
                    int k = IDX(i,j,P.Nx);
                    fprintf(fp,"%g,%g,%g,%g,%g\n",x,y,u[k],v[k],p[k]);
                }
            }
            fclose(fp);
            printf("Saved %s\n", name);
        }
    }

    free(u); free(v); free(p); free(u_star); free(v_star); free(rhs); free(nu_t);
    return 0;
}
