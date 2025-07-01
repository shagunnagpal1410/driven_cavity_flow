#include<iostream>
#include<vector>
#include<map>
#include<cmath>
#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<Eigen/SparseLU>
#include<fstream>
#include<string>
using namespace std;
using namespace Eigen;
typedef Triplet<double> Tri;
struct point{
    double x,y;
    int voxel;
    point(double i, double j) {
        x=i, y=j;
        voxel=0;
    }
    bool operator<(const point& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return voxel < other.voxel;
    }
};
int find_voxel(double x, double y, double x_min, double x_max, int voxels_inrow, double radius) {
    return int((x-x_min)/radius) + int((y)/radius)*voxels_inrow;
}
int find_rownumber(int voxel_number, int voxels_inrow) {
    return voxel_number%voxels_inrow;
}
int find_columnnumber(int voxel_number, int voxels_inrow) {
    return int(voxel_number/voxels_inrow);
}
double gaussian_weight_function(point Ni, point p0, double radius) {
    double distance=(pow(Ni.x-p0.x,2)+pow(Ni.y-p0.y,2))/(radius*radius);
    if (distance<=1) {
        return exp(-6.25*distance);
    }
    else {
        return 0.0;
    }
}
int main() {
    //defining domain
    double L=1.0;
    double B=1.0;
    double dx=0.05;
    double dy=0.05;
    vector<point> previous_domain;
    double x_min=0.0;
    double x_max=1.0;
    double y_min=0.0;
    double y_max=1.0;
    double radius=dx*3;
    int voxels_inrow=int((x_max-x_min)/radius)+1;
    int voxels_incolumn=int(1/radius)+1;
    //initializing some parameters------------------Ending---------------------------------------------------------
    //adding points in the domain-------------------Starting-------------------------------------------------------
    int n=ceil((x_max-x_min)/dx);
    int m=ceil((y_max-y_min)/dy);
    for (int i=0; i<=n; i++) {
        for (int j=0; j<=m; j++) {
            point p1=point(i*dx,j*dy);
            p1.voxel=find_voxel(p1.x,p1.y,x_min,x_max,voxels_inrow,radius);
            previous_domain.push_back(p1);
        }
    }
    int max_voxels=voxels_incolumn*voxels_inrow;
    vector<vector<point>> points_insidevoxel(max_voxels);
    for (int i=0; i<previous_domain.size(); i++) {
        point p0=previous_domain[i];
        points_insidevoxel[p0.voxel].push_back(p0);
    }
    //adding points in the domain-------------------Ending---------------------------------------------------------
    //finding neighbour voxels of each voxel--------Starting-------------------------------------------------------
    vector<vector<int>> neighbour_voxels(max_voxels);
    for(int i=0; i<max_voxels; i++) {
        int find_x=find_rownumber(i, voxels_inrow);
        int find_y=find_columnnumber(i, voxels_inrow);
        for (int diffx=-1; diffx<=1; diffx++) {
            for (int diffy=-1; diffy<=1; diffy++) {
                if (find_x+diffx>=0 && find_x+diffx<=voxels_inrow-1 && diffy+find_y>=0 && find_y+diffy<=voxels_incolumn-1) {
                            neighbour_voxels[i].push_back((find_y+diffy)*voxels_inrow+(find_x+diffx));
                }
            }
        }
    }
    //finding neighbour voxels of each voxel--------Ending--------------------------------------------------------
    //we will be finding neighbours of each point respectively
    map<point,vector<point>> neighbours_ofpoint;
    for (int p=0; p<previous_domain.size(); p++) {
        point p0=previous_domain[p];
        int voxel_number=p0.voxel;
        for (int i=0; i<neighbour_voxels[voxel_number].size(); i++) {
            int neighbour_voxel=neighbour_voxels[voxel_number][i];
            for (int j=0; j<points_insidevoxel[neighbour_voxel].size(); j++) {
                point Ni=points_insidevoxel[neighbour_voxel][j];
                double distance=pow(Ni.x-p0.x,2)+pow(Ni.y-p0.y,2);
                if (distance>0 && distance<pow(radius,2)) {
                    neighbours_ofpoint[p0].push_back(Ni);
                }
            }
        }
    }
    //defining previous velocities
    map<point,double> prev_u;
    map<point,double> prev_v;
    for (int p=0; p<previous_domain.size(); p++) {
        point p0=previous_domain[p];
        if (p0.y==1.0) {
            prev_u[p0]=16*p0.x*p0.x*pow(1-p0.x,2);
            prev_v[p0]=0.0;
        }
        else {
            prev_u[p0]=0.0;
            prev_v[p0]=0.0;
        }
    }
    double T,dt;
    cout<<"Enter the maximum time: ";
    cin>>T;
    cout<<"Enter the time step: ";
    cin>>dt;
    int total_steps=int(T/dt);
    //starting the time loop
    for (int t=1; t<=total_steps; t++) {
        cout<<"time step- "<<t<<" started"<<endl;
        //step-1 the process of calculating v* and u*
        map<point,double> usquare, vsquare, uv, usquarex, vsquarey, uvx, uvy, laplacianu, laplacianv, u_star, v_star;
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            usquare[p0]=pow(prev_u[p0],2.0);
            vsquare[p0]=pow(prev_v[p0],2.0);
            uv[p0]=prev_u[p0]*prev_v[p0];
        }
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            vector<point> neighbours=neighbours_ofpoint[p0];
            int totalneighbours=neighbours.size();
            MatrixXd M(totalneighbours,6);
            MatrixXd W=MatrixXd :: Zero(totalneighbours,totalneighbours);
            VectorXd b1(totalneighbours), b2(totalneighbours), b3(totalneighbours), b4(totalneighbours), b5(totalneighbours);
            for (int i=0; i<totalneighbours; i++) {
                point Ni=neighbours[i];
                double Dx=Ni.x-p0.x, Dy=Ni.y-p0.y;
                M(i,0)=1, M(i,1)=Dx, M(i,2)=Dy, M(i,3)=(Dx*Dx)/2.0, M(i,4)=(Dy*Dy)/2.0, M(i,5)=(Dx*Dy);
                W(i,i)=gaussian_weight_function(Ni,p0,radius);
                b1(i)=prev_u[Ni];
                b2(i)=prev_v[Ni];
                b3(i)=usquare[Ni];
                b4(i)=vsquare[Ni];
                b5(i)=uv[Ni];
            }
            MatrixXd MTWM=M.transpose()*W*M;
            MatrixXd MTW=M.transpose()*W;  
            MatrixXd A=MTWM.ldlt().solve(MTW);
            VectorXd L1(6), L2(6), L3(6);
            L1<<0,0,0,1,1,0;
            L2<<0,1,0,0,0,0;
            L3<<0,0,1,0,0,0;
            VectorXd temp1=A*b1, temp2=A*b2, temp3=A*b3, temp4=A*b4, temp5=A*b5;
            laplacianu[p0]=L1.transpose()*temp1;
            laplacianv[p0]=L1.transpose()*temp2;
            usquarex[p0]=L2.transpose()*temp3;
            vsquarey[p0]=L3.transpose()*temp4;
            uvx[p0]=L2.transpose()*temp5;
            uvy[p0]=L3.transpose()*temp5;
        }
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            double mew=0.01;
            u_star[p0]=prev_u[p0]+(dt*((mew*laplacianu[p0])-(usquarex[p0]+uvy[p0])));
            v_star[p0]=prev_v[p0]+(dt*((mew*laplacianv[p0])-(vsquarey[p0]+uvx[p0])));
        }
        cout<<"We have successfully calculated u* and v*"<<endl;
        //step-2 the process of calculating u*x and v*y
        map<point,double> u_starx, v_stary;
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            vector<point> neighbours=neighbours_ofpoint[p0];
            int totalneighbours=neighbours.size();
            MatrixXd M(totalneighbours,6);
            MatrixXd W=MatrixXd :: Zero(totalneighbours,totalneighbours);
            VectorXd b1(totalneighbours), b2(totalneighbours);
            for (int i=0; i<totalneighbours; i++) {
                point Ni=neighbours[i];
                double Dx=Ni.x-p0.x, Dy=Ni.y-p0.y;
                M(i,0)=1, M(i,1)=Dx, M(i,2)=Dy, M(i,3)=(Dx*Dx)/2.0, M(i,4)=(Dy*Dy)/2.0, M(i,5)=(Dx*Dy);
                W(i,i)=gaussian_weight_function(Ni,p0,radius);
                b1(i)=u_star[Ni];
                b2(i)=v_star[Ni];
            }
            MatrixXd MTWM=M.transpose()*W*M;
            MatrixXd MTW=M.transpose()*W;  
            MatrixXd A=MTWM.ldlt().solve(MTW);
            VectorXd L1(6), L2(6), temp1=A*b1, temp2=A*b2;
            L1<<0,1,0,0,0,0;
            L2<<0,0,1,0,0,0;
            u_starx[p0]=L1.transpose()*temp1;
            v_stary[p0]=L2.transpose()*temp2;
        }
        cout<<"we have successfully calculated v*y and u*x"<<endl;
        //step-3 calculating pressure
        map<point, double> pressure;
        int total_points=previous_domain.size();
        SparseMatrix<double> letsSolve(total_points,total_points);
        VectorXd rhs(total_points);
        vector<Tri> coefficients;
        map<point,int> identity;
        for (int p=0; p<previous_domain.size(); p++) {
            identity[previous_domain[p]]=p;
        }
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            vector<point> neighbours=neighbours_ofpoint[p0];
            int totalneighbours=neighbours.size();
            if ((p0.x==x_min && p0.y==y_min) || (p0.x==x_min && p0.y==y_max) && (p0.x==x_max && p0.y==y_min) || (p0.x==x_max && p0.y==y_max) ) {
                MatrixXd M(totalneighbours+3,6);
                MatrixXd W=MatrixXd :: Zero(totalneighbours+3,totalneighbours+3);
                for (int i=0; i<totalneighbours; i++) {
                    point Ni=neighbours[i];
                    double Dx=Ni.x-p0.x, Dy=Ni.y-p0.y;
                    M(i,0)=1, M(i,1)=Dx, M(i,2)=Dy, M(i,3)=(Dx*Dx)/2.0, M(i,4)=(Dy*Dy)/2.0, M(i,5)=(Dx*Dy);
                    W(i,i)=gaussian_weight_function(Ni,p0,radius);
                }
                M(totalneighbours,0)=0,M(totalneighbours,1)=0,M(totalneighbours,2)=0,M(totalneighbours,3)=1
                ,M(totalneighbours,4)=1,M(totalneighbours,5)=0; 
                W(totalneighbours,totalneighbours)=1;
                M(totalneighbours+1,0)=0,M(totalneighbours+1,1)=1,M(totalneighbours+1,2)=0,M(totalneighbours+1,3)=0
                ,M(totalneighbours+1,4)=0,M(totalneighbours+1,5)=0; 
                W(totalneighbours+1,totalneighbours+1)=1;
                M(totalneighbours+2,0)=0,M(totalneighbours+2,1)=0,M(totalneighbours+2,2)=1,M(totalneighbours+2,3)=0
                ,M(totalneighbours+2,4)=0,M(totalneighbours+2,5)=0; 
                W(totalneighbours+2,totalneighbours+2)=1;
                MatrixXd MTWM=M.transpose()*W*M;
                MatrixXd MTW=M.transpose()*W;  
                MatrixXd A=MTWM.ldlt().solve(MTW);
                double constant=0.0;
                for (int i=0; i<totalneighbours; i++) {
                    point Ni=neighbours[i];
                    coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                }
                double value=(u_starx[p0]+v_stary[p0])/dt;
                constant=constant-(A(0,totalneighbours)*value);
                coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                rhs(identity[p0])=constant;
            }
            else if (p0.x==x_min || p0.x==x_max) {
                MatrixXd M(totalneighbours+2,6);
                MatrixXd W=MatrixXd :: Zero(totalneighbours+2,totalneighbours+2);
                for (int i=0; i<totalneighbours; i++) {
                    point Ni=neighbours[i];
                    double Dx=Ni.x-p0.x, Dy=Ni.y-p0.y;
                    M(i,0)=1, M(i,1)=Dx, M(i,2)=Dy, M(i,3)=(Dx*Dx)/2.0, M(i,4)=(Dy*Dy)/2.0, M(i,5)=(Dx*Dy);
                    W(i,i)=gaussian_weight_function(Ni,p0,radius);
                }
                M(totalneighbours,0)=0,M(totalneighbours,1)=0,M(totalneighbours,2)=0,M(totalneighbours,3)=1
                ,M(totalneighbours,4)=1,M(totalneighbours,5)=0; 
                W(totalneighbours,totalneighbours)=1;
                M(totalneighbours+1,0)=0,M(totalneighbours+1,1)=1,M(totalneighbours+1,2)=0,M(totalneighbours+1,3)=0
                ,M(totalneighbours+1,4)=0,M(totalneighbours+1,5)=0; 
                W(totalneighbours+1,totalneighbours+1)=1;
                MatrixXd MTWM=M.transpose()*W*M;
                MatrixXd MTW=M.transpose()*W;  
                MatrixXd A=MTWM.ldlt().solve(MTW);
                double constant=0.0;
                for (int i=0; i<totalneighbours; i++) {
                    point Ni=neighbours[i];
                    coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                }
                double value=(u_starx[p0]+v_stary[p0])/dt;
                constant=constant-(A(0,totalneighbours)*value);
                coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                rhs(identity[p0])=constant;
            }
            else if (p0.y==y_min || p0.y==y_max) {
                MatrixXd M(totalneighbours+2,6);
                MatrixXd W=MatrixXd :: Zero(totalneighbours+2,totalneighbours+2);
                for (int i=0; i<totalneighbours; i++) {
                    point Ni=neighbours[i];
                    double Dx=Ni.x-p0.x, Dy=Ni.y-p0.y;
                    M(i,0)=1, M(i,1)=Dx, M(i,2)=Dy, M(i,3)=(Dx*Dx)/2.0, M(i,4)=(Dy*Dy)/2.0, M(i,5)=(Dx*Dy);
                    W(i,i)=gaussian_weight_function(Ni,p0,radius);
                }
                M(totalneighbours,0)=0,M(totalneighbours,1)=0,M(totalneighbours,2)=0,M(totalneighbours,3)=1
                ,M(totalneighbours,4)=1,M(totalneighbours,5)=0; 
                W(totalneighbours,totalneighbours)=1;
                M(totalneighbours+1,0)=0,M(totalneighbours+1,1)=0,M(totalneighbours+1,2)=1,M(totalneighbours+1,3)=0
                ,M(totalneighbours+1,4)=0,M(totalneighbours+1,5)=0; 
                W(totalneighbours+1,totalneighbours+1)=1;
                MatrixXd MTWM=M.transpose()*W*M;
                MatrixXd MTW=M.transpose()*W;  
                MatrixXd A=MTWM.ldlt().solve(MTW);
                double constant=0.0;
                for (int i=0; i<totalneighbours; i++) {
                    point Ni=neighbours[i];
                    coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                }
                double value=(u_starx[p0]+v_stary[p0])/dt;
                constant=constant-(A(0,totalneighbours)*value);
                coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                rhs(identity[p0])=constant;
            }
            else {
                MatrixXd M(totalneighbours+1,6);
                MatrixXd W=MatrixXd :: Zero(totalneighbours+1,totalneighbours+1);
                for (int i=0; i<totalneighbours; i++) {
                    point Ni=neighbours[i];
                    double Dx=Ni.x-p0.x, Dy=Ni.y-p0.y;
                    M(i,0)=1, M(i,1)=Dx, M(i,2)=Dy, M(i,3)=(Dx*Dx)/2.0, M(i,4)=(Dy*Dy)/2.0, M(i,5)=(Dx*Dy);
                    W(i,i)=gaussian_weight_function(Ni,p0,radius);
                }
                M(totalneighbours,0)=0,M(totalneighbours,1)=0,M(totalneighbours,2)=0,M(totalneighbours,3)=1
                ,M(totalneighbours,4)=1,M(totalneighbours,5)=0; 
                W(totalneighbours,totalneighbours)=1;
                MatrixXd MTWM=M.transpose()*W*M;
                MatrixXd MTW=M.transpose()*W;  
                MatrixXd A=MTWM.ldlt().solve(MTW);
                double constant=0.0;
                for (int i=0; i<totalneighbours; i++) {
                    point Ni=neighbours[i];
                    coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                }
                double value=(u_starx[p0]+v_stary[p0])/dt;
                constant=constant-(A(0,totalneighbours)*value);
                coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                rhs(identity[p0])=constant;
            }
        }
        SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
        letsSolve.setFromTriplets(coefficients.begin(), coefficients.end());
        solver.compute(letsSolve);
        VectorXd answer = solver.solve(rhs);
        answer=answer.array()-answer.mean();
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            pressure[p0]=answer(identity[p0]);
        }
        cout<<"pressure calculated successfully"<<endl;
        //step-4 calculating px and py
        map<point,double> px,py;
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            vector<point> neighbours=neighbours_ofpoint[p0];
            int totalneighbour=neighbours.size();
            MatrixXd M(totalneighbour,6);
            MatrixXd W=MatrixXd :: Zero(totalneighbour,totalneighbour);
            VectorXd b(totalneighbour);
            for (int i=0; i<totalneighbour; i++) {
                point Ni=neighbours[i];
                double Dx=Ni.x-p0.x, Dy=Ni.y-p0.y;
                M(i,0)=1, M(i,1)=Dx, M(i,2)=Dy, M(i,3)=(Dx*Dx)/2.0, M(i,4)=(Dy*Dy)/2.0, M(i,5)=(Dx*Dy);
                W(i,i)=gaussian_weight_function(Ni,p0,radius);
                b(i)=pressure[Ni];
            }
            MatrixXd MTWM=M.transpose()*W*M;
            MatrixXd MTW=M.transpose()*W;  
            MatrixXd A=MTWM.ldlt().solve(MTW);
            VectorXd L1(6), L2(6), temp=A*b;
            L1<<0,1,0,0,0,0;
            L2<<0,0,1,0,0,0;
            px[p0]=L1.transpose()*temp;
            py[p0]=L2.transpose()*temp;
        }
        cout<<"we have successfully calculated px and py"<<endl;
        //step-5 calculating next velocities
        map<point,double> next_u, next_v;
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            next_u[p0]=u_star[p0]-(dt*px[p0]);
            next_v[p0]=v_star[p0]-(dt*py[p0]);
        }
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            if (p0.y==y_max) {
               next_u[p0]=16*p0.x*p0.x*pow(1-p0.x,2);
               next_v[p0]=0.0;
            }
            else if (p0.y==y_min || p0.x==x_min || p0.x==x_max) {
                next_u[p0]=0.0;
                next_v[p0]=0.0;
            }
            else {
                next_u[p0]=u_star[p0]-dt*px[p0];
                next_v[p0]=v_star[p0]-dt*py[p0];
            }
        }
        prev_u=next_u;
        prev_v=next_v;
        cout<<"time step- "<<t<<" ended successfully"<<endl;
            ofstream fout("Velocity"+to_string(t)+".csv");
            fout<<"X"<<","<<"Y"<<','<<"u"<<","<<"v"<<"\n";
            for (int p=0; p<previous_domain.size(); p++) {
                point p0=previous_domain[p];
                fout<<p0.x<<","<<p0.y<<","<<prev_u[p0]<<","<<prev_v[p0]<<"\n";
        }
        fout.close(); 
    }
    return 0;
}