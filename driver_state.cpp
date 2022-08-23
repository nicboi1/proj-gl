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
    state.image_depth=new float[width*height];
    for(int i=0;i<width*height;i++){
        //state.image_color[i]= 0;
        state.image_color[i]= make_pixel(0,0,0);
        state.image_depth[i]=1;
    }


}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type) {

    std::vector<data_vertex> dv;
    std::vector<data_geometry> dg;
    int i =0;
    while (i < state.num_vertices) {
        data_vertex holder;
        data_geometry holder2;
        holder.data = state.vertex_data + (i * state.floats_per_vertex);
        dv.push_back(holder);
        holder2.data = dv.at(i).data;
        state.vertex_shader(dv.at(i), holder2, state.uniform_data);
        dg.push_back(holder2);
        i++;
    }
    if (type == render_type::triangle) {
        int i = 0;
        while(i <=dg.size()-3) {
            clip_triangle(state, dg.at(i),dg.at(i+1), dg.at(i+2));
            i+=3;
        }

    }
        if(type==render_type::fan){
            int i=0;
            while(i<dg.size()-2){
                clip_triangle(state,dg.at(0),dg.at(i+1),dg.at(2+i));
                i++;
            }

        }
        if(type==render_type::strip){
            int i=0;
            while(i<dg.size()-2){
                clip_triangle(state,dg.at(i),dg.at(i+1),dg.at(i+2));
                i++;
            }
        }
        if(type==render_type::indexed){
            int i=0;
            while(i<dg.size()){
                clip_triangle(state,dg.at(state.index_data[3*1]),dg.at(state.index_data[3*1+1]),dg.at(state.index_data[3*1+2]));
                i++;
            }
        }



}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face) {
    if (face == 6) {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    data_geometry vzed = v0;
    data_geometry vone = v1;
    data_geometry vtwo = v2;
    //for easy readablity
    float v0Z = vzed.gl_Position[2];
    float v0W = vzed.gl_Position[3];
    float v1Z = vone.gl_Position[2];
    float v1W = vone.gl_Position[3];
    float v2Z = vtwo.gl_Position[2];
    float v2W = vtwo.gl_Position[3];

        //one inside
        if (v0Z < -v0W) {
            //out
            if (v1Z >= -v1W) {
                //out
                if (v2Z >= -v2W) {
                   // std::cout<<"Seg2";
                    data_geometry test1[3], test2[3];
                    //lam n points
                    float q = (-v1Z - v2W) / (-v1W + v0W - v1Z + v0Z);
                    vec4 q1 = (vone.gl_Position*q)+(1-q)*vzed.gl_Position;
                    float p = (-v0Z - v0W) / (v2W - v0W + v2Z - v0Z);
                    vec4 p1 = (vtwo.gl_Position*p)+(1-p)*vzed.gl_Position;

                    //set triangle to test
                    test1[0].data = new float[MAX_FLOATS_PER_VERTEX];
                    test1[1] = v1;
                    test1[2] = v2;
                    //interp
                    int i=0;
                    while( i < state.floats_per_vertex) {
                       // std::cout<<"Seg3";
                        interp_type it = state.interp_rules[i];
                        if (it == interp_type::flat) {
                            test1[0].data[i] = vzed.data[i];
                        }else if (it == interp_type::smooth) {
                            test1[0].data[i] = p*vtwo.data[i]+(1-p)*vzed.data[i];
                        }else if (it == interp_type::noperspective) {
                            float o = p*(v2W/(p*v2W+(1-p)*v0W));
                            test1[0].data[i] = o*vtwo.data[i]+(1-o)*vzed.data[i];
                        }
                        i++;
                    }
                    test1[0].gl_Position = p1;
                    //clip with new point [0] and old points
                    clip_triangle(state, test1[0], test1[1], test1[2], face + 1);
                    //set sec triangle
                    test2[0].data = new float[MAX_FLOATS_PER_VERTEX];
                    test2[1] = v1;
                    test2[2] = test1[0];
                    i=0;
                    while (i < state.floats_per_vertex) {
                        //std::cout<<"Seg4";
                        interp_type it = state.interp_rules[i];
                        if (it == interp_type::flat) {
                            test2[0].data[i]=vzed.data[i];
                        } else if (it ==interp_type::smooth) {
                            test2[0].data[i]= q*vzed.data[i]+(1-q)*vone.data[i];
                        }else if (it == interp_type::noperspective) {
                            float o = q*(v0W/(q*v0W+(1-q)*v1W));
                            test2[0].data[i] =o*vzed.data[i]+(1-o)*vone.data[i];
                        }
                        i++;
                    }
                    test2[0].gl_Position = q1;
                    //clip with new new point, new point, and old point
                    clip_triangle(state, test2[0], test2[1], test2[2], face + 1);
                }
            }
        }
        else{
          //  std::cout<<"Seg5";
            //a is out, shuffle
            clip_triangle(state, v1, v2, v0, face+1);
            }
   // }

    }


// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2) {
   // std::cout<<"Seg1rt";
    //std::cout<<"here"<<std::endl;
    vec2 XY1(v0.gl_Position[0]/v0.gl_Position[3],v0.gl_Position[1]/v0.gl_Position[3]);
    vec2 XY2(v1.gl_Position[0]/v1.gl_Position[3],v1.gl_Position[1]/v1.gl_Position[3]);
    vec2 XY3(v2.gl_Position[0]/v2.gl_Position[3],v2.gl_Position[1]/v2.gl_Position[3]);

//each xy after
    vec3 A=vec3(state.image_width/2+XY1[0]*state.image_width/2,state.image_height/2+XY1[1]*state.image_height/2,0.0);
    vec3 B=vec3(state.image_width/2+XY2[0]*state.image_width/2,state.image_height/2+XY2[1]*state.image_height/2,0.0);
    vec3 C=vec3(state.image_width/2+XY3[0]*state.image_width/2,state.image_height/2+XY3[1]*state.image_height/2,0.0);
//box
    int minX = std::min(C[0], (std::min(A[0], B[0])));
    int maxX = std::max(C[0], (std::max(A[0], B[0])));
    int minY = std::min(C[1], (std::min(A[1], B[1])));
    int maxY = std::fmax(C[1], (std::max(A[1], B[1])));
    //must check box is on screen
    if (minX<0){
        minX=0;
    }
    if (minY<0){
        minY=0;
    }
    if (maxX>state.image_width){
        maxX=state.image_width;
    }
    if (maxY>state.image_height){
        maxY=state.image_height;
    }
    /*std::cout<<"maxx:";
    std::cout<<maxX;
    std::cout<<" maxy:";
    std::cout<<maxY;
    std::cout<<" ";*/

    vec3 Z(v0.gl_Position[2]/v0.gl_Position[3],v1.gl_Position[2]/v1.gl_Position[3],v2.gl_Position[2]/v2.gl_Position[3]);

    for (int x = minX; x < maxX; x++) {
        for (int y = minY; y < maxY; y++) {
           /* std::cout<<"minx:";
            std::cout<<minX;
            std::cout<<" miny:";
            std::cout<<minY;
            std::cout<<" ";*/

          /*  std::cout<<"x:";
            std::cout<<x;
            std::cout<<" y:";
            std::cout<<y;
            std::cout<<" ";*/

          //  std::cout<<"Seg2rt";
          // std:: cout<<"here1"<<std::endl;
            vec3 spot=vec3(x,y,0.0);
            //baryc
            //z not the mag
            double total=.5*(cross(B-A,B-C)[2]);
            double alpha=.5*(cross(B-spot,B-C)[2]);
            double beta=.5*(cross(spot-A,spot-C)[2]);
            double gamma=.5*(cross(B-A,B-spot)[2]);

            alpha=alpha/total;
            beta=beta/total;
            gamma=gamma/total;

            double compare=alpha*Z[0]+beta*Z[1]+gamma*Z[2];

            if((alpha>=0) && (beta>=0) && (gamma>=0)) {
                //std::cout<<"Seg3rt";
                data_fragment frag;
                frag.data = new float[MAX_FLOATS_PER_VERTEX];
                data_output outer;
                for (int z = 0; z < state.floats_per_vertex; z++) {
                    // std::cout<<"Seg4rt";
                    //std:: cout<<"here2"<<std::endl;
                    // float smoove;
                    interp_type itype = state.interp_rules[z];
                    if (itype == interp_type::flat) {
                        // std::cout<<"flat";
                        //  std::cout<<"flat ";
                        frag.data[z] = v0.data[z];
                    } else if (itype == interp_type::noperspective) {
                        //  std::cout<<"no";
                        frag.data[z] = v0.data[z]*alpha+v1.data[z]*beta+v2.data[z]*gamma;
                    } else if (itype == interp_type::smooth) {
                        // std::cout<<"smooth";
                        float smoove = (alpha/v0.gl_Position[3]+beta/v1.gl_Position[3]+gamma/v2.gl_Position[3]);
                        float a2, b2, g2;
                        a2 = alpha/(smoove*v0.gl_Position[3]);
                        b2 = beta/(smoove*v1.gl_Position[3]);
                        g2 = gamma/(smoove*v2.gl_Position[3]);
                        frag.data[z] = a2*v0.data[z]+b2*v1.data[z]+g2*v2.data[z];
                    }
                    //std::cout<<"END";
                }
                //  std::cout<<std::endl;
                if(compare<(state.image_depth[x+y*state.image_width])) {
                    state.image_depth[x+y*state.image_width ] = compare;
                    state.fragment_shader(frag, outer, state.uniform_data);
                    state.image_color[x+y*state.image_width] = make_pixel((outer.output_color[0])*255,(outer.output_color[1])*255,(outer.output_color[2])*255);
                }
                //std:: cout<<"here3"<<std::endl;


              //

          }

        }

    }
}






