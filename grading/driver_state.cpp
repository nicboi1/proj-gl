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
        state.image_depth[i]=std::numeric_limits<float>::max();
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
    if (type == render_type::triangle) {
        int i=0;
        while ( i < state.num_vertices) {
            data_vertex holder;
            data_geometry holder2;
            holder.data = state.vertex_data + (i * state.floats_per_vertex);
            dv.push_back(holder);
            holder2.data = dv.at(i).data;
            state.vertex_shader(dv.at(i), holder2, state.uniform_data);
            dg.push_back(holder2);
            i++;
        }
        int j=0;
        while(j <=dg.size()-3) {
                rasterize_triangle(state, dg.at(j),dg.at(j+1), dg.at(j+2));
                j+=3;
        }


    }

}
/*data_geometry inters(data_geometry a,data_geometry b){
    float m=(b.gl_Position[1]-a.gl_Position[1])/(b.gl_Position[0]-a.gl_Position[0]);
    float yinter= a.gl_Position[1]-(m/a.gl_Position[0]);
    if(a.gl_Position[0])
}*/
// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{/*
    // all inside, one outside,two outside, all outside not in view, all outside in view
    array<data_geometry> tresCoord[3];
    tresCoord[0]==v0;
    tresCoord[1]==v1;
    tresCoord[2]==v2;
    bool view[3];
    bool xyView[3][2];
   // int notView=0;
    int fix;
    for(int i=0;i<3;i++){
        if((tresCoord[i].gl_Position[0]>=-1)&&(tresCoord[i].gl_Position[0]<=1) {
            xyView[i][0] = true;
        }else{
            xyView[i][0]=false;
        }
        if((tresCoord[i].gl_Position[1]>=-1)&&(tresCoord[i].gl_Position[1]<=1){
            xyView[i][1]=true;
        }else{
            xyView[i][1]=false;
        }
        if(xyView[i][0] && xyView[i][1]){
            view[i]=true;
        }else{
            view[i]=false;
            fix++;
        }

    }
    if(fix==0){
        //all inside
        face=6;
    }
    else if(fix==1){
        //one outside
        while(view.size())
    }
    else if(fix==2){
        //two outside
    }
    else if(fix=3){
        //all outside
    }*/
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
    //std::cout<<"here"<<std::endl;
    vec2 XY1(v0.gl_Position[0]/v0.gl_Position[3],v0.gl_Position[1]/v0.gl_Position[3]);
    vec2 XY2(v1.gl_Position[0]/v1.gl_Position[3],v1.gl_Position[1]/v1.gl_Position[3]);
    vec2 XY3(v2.gl_Position[0]/v2.gl_Position[3],v2.gl_Position[1]/v2.gl_Position[3]);
    vec3 Z(v0.gl_Position[2]/v0.gl_Position[3],v1.gl_Position[2]/v1.gl_Position[3],v2.gl_Position[2]/v2.gl_Position[3]);


//each xy after
    vec3 A=vec3(state.image_width/2+XY1[0]*state.image_width/2,state.image_height/2+XY1[1]*state.image_height/2,0.0);
    vec3 B=vec3(state.image_width/2+XY2[0]*state.image_width/2,state.image_height/2+XY2[1]*state.image_height/2,0.0);
    vec3 C=vec3(state.image_width/2+XY3[0]*state.image_width/2,state.image_height/2+XY3[1]*state.image_height/2,0.0);
//box
    int minX = fmin(C[0], (fmin(A[0], B[0])));
    int maxX = fmax(C[0], (fmax(A[0], B[0])));
    int minY = fmin(C[1], (fmin(A[1], B[1])));
    int maxY = fmax(C[1], (fmax(A[1], B[1])));

    for (int x = minX; x < maxX; x++) {
        for (int y = minY; y < maxY; y++) {
          // std:: cout<<"here1"<<std::endl;
            vec3 spot=vec3(x,y,0.0);
        //z not the mag
            double total=.5*(cross(B-A,B-C)[2]);
            double alpha=.5*(cross(B-spot,B-C)[2]);
            double beta=.5*(cross(spot-A,spot-C)[2]);
            double gamma=.5*(cross(B-A,B-spot)[2]);

            alpha=alpha/total;
            beta=beta/total;
            gamma=gamma/total;

            double compare=alpha*Z[0]+beta*Z[1]+gamma*Z[2];

            if((alpha>=0) && (beta>=0) && (gamma>=0)){
                data_fragment frag;
                frag.data=new float[MAX_FLOATS_PER_VERTEX];
                data_output outer;

                // A=v B=xyz C=n

              for (int k = 0; k < state.floats_per_vertex; k++) {
                  //std:: cout<<"here2"<<std::endl;

                  interp_type itype=state.interp_rules[k];
                  if(itype==interp_type::flat) {
                    //  std::cout<<"flat ";
                      frag.data[k] = v0.data[k];
                  }
                  else if(itype==interp_type::noperspective){
                    //  std::cout<<"no";
                      frag.data[k] = v0.data[k] * alpha + v1.data[k] * beta + v2.data[k] * gamma;
                      }
                  else if(itype==interp_type::smooth){

                  }
              }
              //  std::cout<<std::endl;

              state.fragment_shader(frag, outer, state.uniform_data);

              if(compare<(state.image_depth[y*state.image_width+x])) {
                  state.image_depth[y*state.image_width+x]=compare;
                   //std:: cout<<"here3"<<std::endl;
                  state.image_color[y*state.image_width+x] = make_pixel((outer.output_color[0]) * 255,
                                                                            (outer.output_color[1]) * 255,
                                                                            (outer.output_color[2]) * 255);
                }
          }

        }

    }
}






