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
        for (int j = 0; j <= all.size()-3; j+=3) {

                rasterize_triangle(state, all.at(j),all.at(j+1), all.at(j+2));

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
    vec3 VT1(v0.gl_Position[0]/v0.gl_Position[3],v0.gl_Position[1]/v0.gl_Position[3],v0.gl_Position[2]/v0.gl_Position[3]);
    vec3 VT2(v1.gl_Position[0]/v1.gl_Position[3],v1.gl_Position[1]/v1.gl_Position[3],v1.gl_Position[2]/v1.gl_Position[3]);
    vec3 VT3(v2.gl_Position[0]/v2.gl_Position[3],v2.gl_Position[1]/v2.gl_Position[3],v2.gl_Position[2]/v2.gl_Position[3]);

    float Ax=(state.image_width/2)+VT1[0]*(state.image_width/2);
    float Ay=(state.image_height/2)+VT1[1]*(state.image_height/2);
    float Bx=(state.image_width/2)+VT2[0]*(state.image_width/2);
    float By=(state.image_height/2)+VT2[1]*(state.image_height/2);
    float Cx=(state.image_width/2)+VT3[0]*(state.image_width/2);
    float Cy=(state.image_height/2)+VT3[1]*(state.image_height/2);
    //float zbuf=std::numeric_limits<max>;

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

        if((alpha>0) && (beta>0)  && (gamma>0) && (alpha+beta+gamma<=1+1e-4)){
            data_fragment frag;
            data_output outer;
            //inside triangle
           // A=v B=xyz C=n
            //float z=alpha*v0.gl_Position[2]+beta*v1.gl_Position[2]*gamma*v2.gl_Position[2];




                state.fragment_shader(frag,outer,state.uniform_data);
                state.image_color[y*state.image_width+x]= make_pixel( (outer.output_color[0])*255, (outer.output_color[1])*255, (outer.output_color[2])*255);


            }
        }
    }



    }


