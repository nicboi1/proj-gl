#include "driver_state.h"
#include <cstring>
#include <vector>
driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=new pixel[width*height];
    for(int i=0;i<width*height;i++){
        //state.image_color[i]= 0;
        state.image_color[i]= make_pixel(0,0,0);
    }
    state.image_depth=new float[width*height];
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type) {

    std::vector <data_geometry> all;
    std::vector <data_vertex> aye;
    if (type == render_type::triangle) {

        int j=0;
        while(j<state.num_vertices){
            data_geometry trial;
            data_vertex tryme;
            tryme.data=state.vertex_data+(j*state.floats_per_vertex);
            aye.push_back(tryme);
            state.vertex_shader(aye.at(j),trial,state.uniform_data);
            all.push_back(trial);
            j++;
        }
        for (int j = 0; j < all.size(); j++) {
            if (j % 3 == 0) {
                rasterize_triangle(state, all[j], all[j+1], all[j+2]);
            }
        }


    }

}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,v0,v1,v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2) {

    float Ax=(state.image_width/2)+v0.gl_Position[0]*(state.image_width/2);
    float Ay=(state.image_height/2)+v0.gl_Position[1]*(state.image_height/2);
    float Bx=(state.image_width/2)+v1.gl_Position[0]*(state.image_width/2);
    float By=(state.image_height/2)+v1.gl_Position[1]*(state.image_height/2);
    float Cx=(state.image_width/2)+v2.gl_Position[0]*(state.image_width/2);
    float Cy=(state.image_height/2)+v2.gl_Position[1]*(state.image_height/2);


    int minX = fmin(Cx, (fmin(Ax, Bx)));
    int maxX = fmax(Cx, (fmax(Ax, Bx)));
    int minY = fmin(Cy, (fmin(Ay, By)));
    int maxY = fmax(Cy, (fmax(Ay, By)));



    for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {

        //v0=a v1=b v2=c ij=p


     float total= (Ax*By+Bx*Cy+Cx*Ay-Ay*Bx-By*Cx-Cy*Ax)/2;
     float alpha= (x*By+Bx*Cy+Cx*y-y*Bx-By*Cx-Cy*x)/2;
     float beta= (Ax*y+x*Cy+Cx*Ay-Ay*x-y*Cx-Cy*Ax)/2;
     float gamma= (Ax*By+Bx*y+x*Ay-Ay*Bx-By*x-y*Ax)/2;

        alpha=alpha/total;
        beta=beta/total;
        gamma=gamma/total;

        if((alpha>0) && (beta>0)  && (gamma>0) && (alpha+beta+gamma<=1)){
            //inside triangle
          /*  float c0[3]={Ax.data[3],v0.data[4],v0.data[5]};
            float c1[3]={v1.data[3],v1.data[4],v1.data[5]};
            float c2[3]={v2.data[3],v2.data[4],v2.data[5]};

            vec3 c00={alpha*c0[0],beta*c0[1],gamma*c0[2]};
            vec3 c01={alpha*c1[0],beta*c1[1],gamma*c1[2]};
            vec3 c02={alpha*c2[0],beta*c2[1],gamma*c2[2]};

            vec3 c00={alpha*c0[0],alpha*c0[1],alpha*c0[2]};
            vec3 c01={beta*c1[0],beta*c1[1],beta*c1[2]};
            vec3 c02={gamma*c2[0],gamma*c2[1],gamma*c2[2]};

            float color[3];
            color[0]=c00[0]+c01[0]+c02[0];
            color[1]=c00[1]+c01[1]+c02[1];
            color[2]=c00[2]+c01[2]+c02[2];
            */
            state.image_color[y*state.image_width+x]= make_pixel(255,255,255);
            }
        }
    }



    }


