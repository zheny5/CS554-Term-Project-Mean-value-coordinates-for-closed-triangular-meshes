/*

Functions for learnply

Eugene Zhang, 2005
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "glut.h"
#include <string.h>
#include <fstream>
#include "icVector.H"
#include "icMatrix.H"
#include "learnply.h"
#include "trackball.h"
#include "tmatrix.h"
#include  <vector>
Polyhedron* poly;
Polyhedron* poly2;
std::vector<icVector3> poly2_colors;//to set the color for control mesh

const int win_width  = 1024;
const int win_height = 1024;

bool display_wireframe = false;
bool display_control_mesh_frame = false;
bool display_control_mesh = false;

int  display_mode = 0;
int view_mode = 0;   // 0 = othogonal, 1=perspective
int ACSIZE = 1;      // for antialiasing

int  last_x = 0, last_y = 0;
double radius_factor = 0.9;

float s_old, t_old;
float rotmat[4][4];
static Quaternion rvec;

int mouse_mode = -2;   // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right

struct jitter_struct { double x; double y; };
jitter_struct ji1[1] = { {0.0, 0.0} };
jitter_struct ji16[16] = {{0.125, 0.125}, {0.375, 0.125}, {0.625, 0.125}, {0.875, 0.125},
						  {0.125, 0.375}, {0.375, 0.375}, {0.625, 0.375}, {0.875, 0.375},
						  {0.125, 0.625}, {0.375, 0.625}, {0.625, 0.625}, {0.875, 0.625},
						  {0.125, 0.875}, {0.375, 0.875}, {0.625, 0.875}, {0.875, 0.875} };
float L_checker = 0.5;//for the checkerboard

void init(void);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void display_shape(GLenum mode, Polyhedron* poly);

/******************************************************************************
Main program.
******************************************************************************/
int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		FILE* this_file = fopen("../tempmodels/dog.ply", "r");
		if (this_file == NULL) { return -1; }
		poly = new Polyhedron(this_file);
		fclose(this_file);
		//load control mesh
		FILE* this_file2 = fopen("../tempmodels/dog_control_mesh.ply", "r");
		if (this_file2 == NULL) { return -1; }
		poly2 = new Polyhedron(this_file2);
		fclose(this_file2);
	}
	else
	{ 
		FILE* this_file = fopen(argv[1], "r");
		if (this_file == NULL) { return -1; }
		poly = new Polyhedron(this_file);
		fclose(this_file);
		//load control mesh
		FILE* this_file2 = fopen(argv[2], "r");
		if (this_file2 == NULL) { return -1; }
		poly2 = new Polyhedron(this_file2);
		fclose(this_file2);
	}
	poly->initialize(); // initialize everything
	poly2->initialize(); // initialize everything
	poly2->generateColorControlMesh();//generate color for control mesh
	std::vector<std::vector<Weights>> meanValueCoordinates = poly->getMeanValueCoordinates(*poly2);
	poly->InterpolateColor(*poly2, *poly, meanValueCoordinates);//Interpolate color for poly(model mesh) based on poly2(control mesh)
	//poly->InterpolatePos(*poly2, *poly, meanValueCoordinates);//Interpolate position for deformation on poly(model mesh) based on poly2(control mesh)

	mat_ident(rotmat);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("Geometric Modeling");
	init();
	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMainLoop();
	poly->finalize();  // finalize everything
	poly2->finalize();// finalize everything

	return 0;    /* ANSI C requires main to return int. */
}

int Parity(int n) {
	return n % 2 == 0 ? 1 : 0;
}
void init(void) {
	/* select clearing color */
	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	// may need it
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glEnable(GL_NORMALIZE);
	if (poly->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/
void keyboard(unsigned char key, int x, int y) {
	int i;

	/* set escape key to exit */
	switch (key) {
	case 27:
		poly->finalize();  // finalize_everything
		poly2->finalize();// finalize_everything
		exit(0);
		break;

	case '1':
		display_wireframe = !display_wireframe;
		display();
		break;

	case '2':
		display_mode = 0;
		display();
		break;
	case '3':
		display_mode = 3;//checker board, either subdivision or not
		display();
		break;
	case '4':
		display_mode = 4;//display the frame for control mesh
		display_control_mesh_frame = !display_control_mesh_frame;
		display();
		break;
	case '5'://display colored mesh
		display_mode = 5;
		display();
		break;
	case '6'://display colored control mesh
		display_mode = 6;
		display_control_mesh = !display_control_mesh;
		display();
		break;

	case 'x':
		switch (ACSIZE) 
		{
		case 1:
			ACSIZE = 16;
			break;
		default:
			ACSIZE = 1;
			break;
		}
		fprintf(stderr, "ACSIZE=%d\n", ACSIZE);
		display();
		break;

	case '|':
	{
		FILE* this_file = fopen("rotmat.txt", "w");
		for (i = 0; i < 4; i++)
			fprintf(this_file, "%f %f %f %f\n", rotmat[i][0], rotmat[i][1], rotmat[i][2], rotmat[i][3]);
		fclose(this_file);
	}
	break;

	case '^':
	{
		FILE* this_file = fopen("rotmat.txt", "r");
		if (this_file == NULL) { break; }
		for (i = 0; i < 4; i++)
			fscanf(this_file, "%f %f %f %f ", (&rotmat[i][0]), (&rotmat[i][1]), (&rotmat[i][2]), (&rotmat[i][3]));
		fclose(this_file);
		display();
	}
	break;

	}
}

void multmatrix(const Matrix m)
{
	int i, j, index = 0;

	GLfloat mat[16];

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			mat[index++] = m[i][j];

	glMultMatrixf(mat);
}

void set_view(GLenum mode, Polyhedron* poly)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = { 0.3f, 0.3f, 0.3f, 1.0f };
	GLfloat light_diffuse0[] = { 0.7f, 0.7f, 0.7f, 1.0f };
	GLfloat light_specular0[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat light_ambient1[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat light_diffuse1[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	GLfloat light_specular1[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat light_ambient2[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat light_diffuse2[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat light_specular2[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat light_position[] = { 0.0f, 0.0f, 0.0f, 1.0f };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);

	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor, radius_factor, -radius_factor, radius_factor, 0.0, 40.0);
	else
		gluPerspective(45.0, 1.0, 0.1, 40.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5f;
	light_position[1] = 0.0f;
	light_position[2] = 0.0f;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1f;
	light_position[1] = 0.0f;
	light_position[2] = 0.0f;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

void set_scene(GLenum mode, Polyhedron* poly)
{
	glTranslatef(0.0, 0.0, -3.0);
	multmatrix(rotmat);

	glScalef(1.0 / poly->radius, 1.0 / poly->radius, 1.0 / poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

void motion(int x, int y) {
	float r[4];
	float xsize, ysize, s, t;

	switch (mouse_mode) {
	case -1:

		xsize = (float)win_width;
		ysize = (float)win_height;

		s = (2.0 * x - win_width) / win_width;
		t = (2.0 * (win_height - y) - win_height) / win_height;

		if ((s == s_old) && (t == t_old))
			return;

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

int processHits(GLint hits, GLuint buffer[])
{
	GLuint names, * ptr;
	double smallest_depth = 1.0e+20, current_depth;
	int seed_id = -1;
	unsigned char need_to_update;

	printf("hits = %d\n", hits);
	ptr = (GLuint*)buffer;
	for (int i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double)*ptr / 0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (unsigned int j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	printf("triangle id = %d\n", seed_id);
	return seed_id;
}

void mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON)
	{
		switch (mouse_mode) {
		case -2:  // no action
			if (state == GLUT_DOWN) {
				float xsize = (float)win_width;
				float ysize = (float)win_height;

				float s = (2.0 * x - win_width) / win_width;
				float t = (2.0 * (win_height - y) - win_height) / win_height;

				s_old = s;
				t_old = t;

				mouse_mode = -1;  // down
				mouse_button = button;
				last_x = x;
				last_y = y;
			}
			break;

		default:
			if (state == GLUT_UP) {
				button = -1;
				mouse_mode = -2;
			}
			break;
		}
	}
	else if (button == GLUT_MIDDLE_BUTTON)
	{
		if (state == GLUT_DOWN)
		{  // build up the selection feedback mode

			GLuint selectBuf[win_width];
			GLint hits;
			GLint viewport[4];

			glGetIntegerv(GL_VIEWPORT, viewport);

			glSelectBuffer(win_width, selectBuf);
			(void)glRenderMode(GL_SELECT);

			glInitNames();
			glPushName(0);

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
			glLoadIdentity();
			/*  create 5x5 pixel picking region near cursor location */
			gluPickMatrix((GLdouble)x, (GLdouble)(viewport[3] - y),
				1.0, 1.0, viewport);

			set_view(GL_SELECT, poly);
			glPushMatrix();
			set_scene(GL_SELECT, poly);
			display_shape(GL_SELECT, poly);
			glPopMatrix();
			glFlush();

			hits = glRenderMode(GL_RENDER);
			poly->seed = processHits(hits, selectBuf);
			display();
		}
	}
}

void display_object()
{
	Polyhedron* the_patch = poly;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	for (int i = 0; i < poly->ntris(); i++)
	{
		Triangle* temp_t = poly->tlist[i];
		glBegin(GL_POLYGON);
		GLfloat mat_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };

		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

		glColor3f(1.0, 1.0, 1.0);
		glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
		for (int j = 0; j < 3; j++)
		{
			Vertex* temp_v = temp_t->verts[j];
			glVertex3d(temp_v->pos.x, temp_v->pos.y, temp_v->pos.z);
		}
		glEnd();
	}
}

void display_shape(GLenum mode, Polyhedron* this_poly)
{
	GLfloat mat_diffuse[4];

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	//render faces
	for (int i = 0; i < this_poly->ntris(); i++)
	{
		if (mode == GL_SELECT) { glLoadName(i + 1); }
		Triangle* temp_t = this_poly->tlist[i];

		if (i == this_poly->seed)
		{
			mat_diffuse[0] = 0.0;
			mat_diffuse[1] = 0.0;
			mat_diffuse[2] = 1.0;
			mat_diffuse[3] = 1.0;
		}
		else {
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			mat_diffuse[3] = 1.0;
		}
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);


		glBegin(GL_POLYGON);
		for (int j = 0; j < 3; j++)
		{
			Vertex* temp_v = temp_t->verts[j];
			glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
			if (i == this_poly->seed)
			{
				glColor3f(0.0, 0.0, 1.0);
			}
			else
			{
				switch (display_mode)
				{
				case 1://
					glColor3f(1.0, 0.0, 0.0);
					break;
				case 3:
					glColor3d(Parity(floor(temp_v->pos.x / L_checker)), Parity(floor(temp_v->pos.y / L_checker)), Parity(floor(temp_v->pos.z / L_checker)));
					glVertex3d(temp_v->pos.x, temp_v->pos.y, temp_v->pos.z);
					break;
				case 5:
					glColor3d(temp_v->color.x, temp_v->color.y, temp_v->color.z);
					glVertex3d(temp_v->pos.x, temp_v->pos.y, temp_v->pos.z);
					break;
				default:
					glColor3f(1.0, 1.0, 0.0);
					break;
				}
			}
			glVertex3d(temp_v->pos.x, temp_v->pos.y, temp_v->pos.z);
		}
		glEnd();
	}

	if (display_wireframe && mode != GL_SELECT)
	{
		glDisable(GL_LIGHTING);
		glLineWidth(2);
		glColor3f(0.0, 0.0, 0.0);
		for (int i = 0; i < this_poly->nedges(); i++)
		{
			Edge* temp_e = this_poly->elist[i];
			glBegin(GL_LINES);
			glVertex3d(temp_e->verts[0]->pos.x, temp_e->verts[0]->pos.y, temp_e->verts[0]->pos.z);
			glVertex3d(temp_e->verts[1]->pos.x, temp_e->verts[1]->pos.y, temp_e->verts[1]->pos.z);
			glEnd();
		}
	}
	
}

void display_control_frame(Polyhedron* this_poly) 
{

	if(display_control_mesh)
	{
		display_shape(GL_RENDER, this_poly);//display control mesh
	
	}
	if (display_control_mesh_frame)
	{
		glDisable(GL_LIGHTING);
		glLineWidth(2);
		glColor3f(0.0, 0.0, 0.0);
		for (int i = 0; i < this_poly->nedges(); i++)
		{
			Edge* temp_e = this_poly->elist[i];
			glBegin(GL_LINES);
			glVertex3d(temp_e->verts[0]->pos.x, temp_e->verts[0]->pos.y, temp_e->verts[0]->pos.z);
			glVertex3d(temp_e->verts[1]->pos.x, temp_e->verts[1]->pos.y, temp_e->verts[1]->pos.z);
			glEnd();
		}
	}

}

void display(void)
{
	GLint viewport[4];
	int jitter;

	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glGetIntegerv(GL_VIEWPORT, viewport);

	glClear(GL_ACCUM_BUFFER_BIT);
	for (jitter = 0; jitter < ACSIZE; jitter++) 
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		set_view(GL_RENDER, poly);
		glPushMatrix();
		switch (ACSIZE) {
		case 1:
			glTranslatef(ji1[jitter].x * 2.0 / viewport[2], ji1[jitter].y * 2.0 / viewport[3], 0.0);
			break;

		case 16:
			glTranslatef(ji16[jitter].x * 2.0 / viewport[2], ji16[jitter].y * 2.0 / viewport[3], 0.0);
			break;

		default:
			glTranslatef(ji1[jitter].x * 2.0 / viewport[2], ji1[jitter].y * 2.0 / viewport[3], 0.0);
			break;
		}
		set_scene(GL_RENDER, poly);
		display_shape(GL_RENDER, poly);//display control mesh
		display_control_frame(poly2);
		glPopMatrix();
		glAccum(GL_ACCUM, 1.0 / ACSIZE);
	}
	glAccum(GL_RETURN, 1.0);
	glFlush();
	glutSwapBuffers();
	glFinish();
}