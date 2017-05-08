/*Matthew Chiborak
250748631
mchibora@uwo.ca
CS3388 Assignment 2
October 28, 2016*/

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else 
#include <GL/glut.h>
#endif

#define ROTATION_ANGLE 2 //Rotate the plane of vertices this value degrees with each transformation
#define PI 3.14159265

static GLfloat LIGHT_ROTATION = 1; //How much the light rotates every incremement
static GLfloat CURRENT_LIGHT_ANGLE = 0; //The current rotation of the light
static GLfloat startLightPosition[] = { 10.0, -30.0, 0.0, 1.0 }; //Starting position
static GLfloat currentLightPosition[] = { 10.0, -30.0, 0.0, 1.0 }; //Lights current position
static GLdouble cameraPosition[] = { 0, 45, -5, 0, 0, 15, 0, 0, -1 }; //Current camera position
static bool viewingMode = false; //Used to deside to show either wireframe or model
static GLfloat viewingProjection[] = { 60.0, 4.0 / 3.0, 1, 80 };

//Set up the different types of light. Note to self. Have all 3 values be the same to create different intensity white light. 
//Also note to self. Amibient + Diffuse should equal 1.0
GLfloat ambientLight[] = { 0.2, 0.2, 0.2, 1.0 };
GLfloat diffuseLight[] = { 0.8, 0.8, 0.8, 1.0 };
GLfloat specularLight[] = { 1.0, 1.0, 1.0, 1.0 };

//Set up the materials
GLfloat vaseShine = 10.0; //0-128
GLfloat vaseColour[] = { 1.0, 0.0, 0.5, 1.0 };
GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat red[] = { 1.0, 0.0, 0.0, 1.0 };
GLfloat green[] = { 0.0, 1.0, 0.0, 1.0 };
GLfloat blue[] = { 0.0, 1.0, 1.0, 1.0 };

//Class for holding onto the xyz coords of the vertices
class Vertex
{
	public:
		float x, y, z, homogeneous;

		Vertex(float x, float y, float z, float homogeneous)
		{
			this->x = x;
			this->y = y;
			this->z = z;
			this->homogeneous = homogeneous;
		}
};

//Class for storing which vertices belong to what polygon face
class Face
{
	public:
		Vertex* v1;
		Vertex* v2;
		Vertex* v3; //Remember these go in counterclockwise
		float normal[3]; //The face's normal. Used for lighting

		Face(Vertex* v1, Vertex* v2, Vertex* v3)
		{
			this->v1 = v1;
			this->v2 = v2;
			this->v3 = v3;

			//Find the normal
			//Find the polygon's normal
			//n = axb / ||axb||
			float aCrossB[3];
			aCrossB[0] = (v2->y - v1->y) * (v3->z - v1->z) - (v2->z - v1->z) * (v3->y - v1->y);
			aCrossB[1] = (v2->z - v1->z) * (v3->x - v1->x) - (v2->x - v1->x) * (v3->z - v1->z);
			aCrossB[2] = (v2->x - v1->x) * (v3->y - v1->y) - (v2->y - v1->y) * (v3->x - v1->x);
			float magnitude = sqrtf(pow(aCrossB[0], 2) + pow(aCrossB[1], 2) + pow(aCrossB[2], 2));

			normal[0] = aCrossB[0] / magnitude;
			normal[1] = aCrossB[1] / magnitude;
			normal[2] = aCrossB[2] / magnitude;
		}
};

//Class to store each rotation of the file read profile
class Profile
{
	public:
		std::vector<Vertex*> vertexHolder;

		void addVertex(Vertex* newVertex)
		{
			vertexHolder.push_back(newVertex);
		}
};

//Class to store all the information about the edges made from triangulation
class Edge
{
	public:
		Vertex* point1;
		Vertex* point2;

		Edge(Vertex* point1, Vertex* point2)
		{
			this->point1 = point1;
			this->point2 = point2;
		}
};

//Vectors to store all the mesh information
std::vector<Profile*> profileHolder;
std::vector<Face*> faceHolder;
std::vector<Edge*> edgeHolder;

//Take the read profile and rotate it and create a new profile. Repeats until full rotation is completed.
void generateRotationProfiles()
{
	//Start at 1st incremement because 1st profile is already done
	int currentAngle = ROTATION_ANGLE;
	int currentProfile = 1;

	while (currentAngle < 360)
	{
		//Create the new profile
		profileHolder.push_back(new Profile());

		//Rotate each point from the original profile around the z-axis and add them to the newly created profile
		//Rotate around the z-axis. x' = xcos - ysin. y' = xsin + ycos. z' = z
		//Iterate over all the vertices in the profile
		for (int i = 0; i < profileHolder.at(0)->vertexHolder.size(); i++)
		{
			float tempx = (profileHolder.at(0)->vertexHolder.at(i)->x) * cos(PI * currentAngle / 180.0) - (profileHolder.at(0)->vertexHolder.at(i)->y) * sin(PI * currentAngle / 180.0);
			float tempy = (profileHolder.at(0)->vertexHolder.at(i)->x) * sin(PI * currentAngle / 180.0) + (profileHolder.at(0)->vertexHolder.at(i)->y) * cos(PI * currentAngle / 180.0);
			float tempz = profileHolder.at(0)->vertexHolder.at(i)->z;
			float temph = profileHolder.at(0)->vertexHolder.at(i)->homogeneous;
			profileHolder.at(currentProfile)->addVertex(new Vertex(tempx, tempy, tempz, temph));
		}
		
		//Set up the next pass
		currentProfile++;
		currentAngle += ROTATION_ANGLE;
	}
}

//Create the wire mesh. Fills the edgeholder with edges between the vertices
void triangulatePoints()
{
	int currentProfile = 0;

	//Check if theres more than one profile to use
	if (profileHolder.size() < 2)
	{
		return;
	}

	//Iterate over all the profiles
	for (int i = 0; i < profileHolder.size(); i++)
	{
		int secondProfileId = i + 1;
		//Account for if the 2nd profile is actually the original and we've looped all the way around
		if (secondProfileId >= profileHolder.size())
		{
			secondProfileId = 0;
		}

		//Minus 1 because the last set of edges will be found with the 2nd last iteration
		for (int j = 0; j < profileHolder.at(i)->vertexHolder.size() - 1; j++)
		{
			//Create the edges. Top, left, bottom, right, center
			edgeHolder.push_back(new Edge(profileHolder.at(i)->vertexHolder.at(j), profileHolder.at(secondProfileId)->vertexHolder.at(j)));
			edgeHolder.push_back(new Edge(profileHolder.at(i)->vertexHolder.at(j), profileHolder.at(i)->vertexHolder.at(j+1)));
			edgeHolder.push_back(new Edge(profileHolder.at(i)->vertexHolder.at(j+1), profileHolder.at(secondProfileId)->vertexHolder.at(j+1)));
			edgeHolder.push_back(new Edge(profileHolder.at(secondProfileId)->vertexHolder.at(j), profileHolder.at(secondProfileId)->vertexHolder.at(j+1)));
			edgeHolder.push_back(new Edge(profileHolder.at(i)->vertexHolder.at(j), profileHolder.at(secondProfileId)->vertexHolder.at(j+1)));

			//Create the face. Remember, counterclockwise
			faceHolder.push_back(new Face(profileHolder.at(i)->vertexHolder.at(j), profileHolder.at(i)->vertexHolder.at(j + 1), profileHolder.at(secondProfileId)->vertexHolder.at(j)));
			faceHolder.push_back(new Face(profileHolder.at(secondProfileId)->vertexHolder.at(j + 1), profileHolder.at(secondProfileId)->vertexHolder.at(j), profileHolder.at(i)->vertexHolder.at(j + 1)));
		}
	}
}

//Function for reading the point file. Returns false if cannot open the file. Takes a string that is the name of the text file
bool readVertices(std::string fileName)
{
	std::string line;
	std::ifstream vertexFile(fileName);
	std::string tempAtt[4];

	//Create the 1st profile
	profileHolder.push_back(new Profile());

	if (vertexFile.is_open())
	{
		while (getline(vertexFile, line))
		{
			//To fix issue where homogeneous coord is being ignored
			line = line + " ";

			std::stringstream ss(line);
			std::string temp = "";
			char i;
			int currentAtt = 0;

			while (ss >> i)
			{
				temp += i;
				char nextChar = ss.peek();

				if (nextChar == ' ')
				{
					ss.ignore();
					tempAtt[currentAtt++] = temp;

					//Reset the holder of the temporary value
					temp = "";
				}
			}

			//Reached the end of a line. Create and store the vertex
			profileHolder.at(0)->addVertex(new Vertex(std::stof(tempAtt[0]), std::stof(tempAtt[1]), std::stof(tempAtt[2]), std::stof(tempAtt[3])));
		}
	}
	else
	{
		//File couldn't be read
		return false;
	}

	//File successfully read and stored vertices
	return true;
}

//Spin the light around the vase and also change the camera position if a change was input
void spinLight()
{
	glPushMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(viewingProjection[0], viewingProjection[1], viewingProjection[2], viewingProjection[3]);

	//Apply the light effects
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);

	//Spin the light
	CURRENT_LIGHT_ANGLE = (CURRENT_LIGHT_ANGLE + LIGHT_ROTATION);

	if (CURRENT_LIGHT_ANGLE >= 360)
	{
		CURRENT_LIGHT_ANGLE = 0;
	}

	//Find the light's new position based on the new angle
	float tempx = (startLightPosition[0]) * cos(PI * CURRENT_LIGHT_ANGLE / 180.0) - (startLightPosition[1]) * sin(PI * CURRENT_LIGHT_ANGLE / 180.0);
	float tempy = (startLightPosition[0]) * sin(PI * CURRENT_LIGHT_ANGLE / 180.0) + (startLightPosition[1]) * cos(PI * CURRENT_LIGHT_ANGLE / 180.0);
	currentLightPosition[0] = tempx;
	currentLightPosition[1] = tempy; 

	//Set the light to its new position
	glLightfv(GL_LIGHT0, GL_POSITION, currentLightPosition);
	glPopMatrix();

	glutPostRedisplay(); //Call display

	//Change the camera position if there was a change
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraPosition[0], cameraPosition[1], cameraPosition[2], cameraPosition[3], cameraPosition[4], cameraPosition[5], cameraPosition[6], cameraPosition[7], cameraPosition[8]);
}


// Clears the window and draw the vase
void display() {

	//Clear the screen
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	//Apply the material
	glMaterialfv(GL_FRONT, GL_AMBIENT, vaseColour);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, vaseColour);
	glMaterialfv(GL_FRONT, GL_SPECULAR, white);
	glMaterialf(GL_FRONT, GL_SHININESS, vaseShine); 

	//If in wireframe mode, draw the wireframe. If not, draw the faces
	if (viewingMode)
	{
		//Draw all the edges to create a wireframe
		for (int i = 0; i < edgeHolder.size(); i++)
		{
			glBegin(GL_LINES);
			glVertex3f(edgeHolder.at(i)->point1->x, edgeHolder.at(i)->point1->y, edgeHolder.at(i)->point1->z); // origin of the line
			glVertex3f(edgeHolder.at(i)->point2->x, edgeHolder.at(i)->point2->y, edgeHolder.at(i)->point2->z); // ending point of the line
			glEnd();
		}
	}
	else
	{
		//Draw all the faces to create the vase
		for (int i = 0; i < faceHolder.size(); i++)
		{
			glBegin(GL_TRIANGLES);
			glNormal3f(faceHolder.at(i)->normal[0], faceHolder.at(i)->normal[1], faceHolder.at(i)->normal[2]);
			glVertex3f(faceHolder.at(i)->v1->x, faceHolder.at(i)->v1->y, faceHolder.at(i)->v1->z);
			glVertex3f(faceHolder.at(i)->v2->x, faceHolder.at(i)->v2->y, faceHolder.at(i)->v2->z);
			glVertex3f(faceHolder.at(i)->v3->x, faceHolder.at(i)->v3->y, faceHolder.at(i)->v3->z);
			glEnd();
		}
	}

	//Draw the axies and rotate them based on the light's rotation
	glPushMatrix();

	glRotatef(CURRENT_LIGHT_ANGLE, 0, 0, 1);

	glBegin(GL_LINES);

	glMaterialfv(GL_FRONT, GL_AMBIENT, red);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
	glVertex3f(0, 0, 0); 
	glVertex3f(100, 0, 0);

	glMaterialfv(GL_FRONT, GL_AMBIENT, green);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 100, 0);

	glMaterialfv(GL_FRONT, GL_AMBIENT, blue);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, 100);

	glEnd();

	glPopMatrix();

	glutSwapBuffers();
}

//The screen was resized
void reshape(int w, int h)
{
	
}

//Release the nemory used to create the vertices, profiles, egdes, and faces
void releaseMemory()
{
	//Release memory
	for (int i = 0; i < profileHolder.size(); i++)
	{
		if (profileHolder.at(i) != NULL)
		{
			for (int j = 0; j < profileHolder.at(i)->vertexHolder.size(); j++)
			{
				if (profileHolder.at(i)->vertexHolder.at(j) != NULL)
				{
					delete profileHolder.at(i)->vertexHolder.at(j);
				}
			}
			delete profileHolder.at(i);
		}
	}
	for (int i = 0; i < edgeHolder.size(); i++)
	{
		if (edgeHolder.at(i) != NULL)
		{
			delete edgeHolder.at(i);
		}
	}
	for (int i = 0; i < faceHolder.size(); i++)
	{
		if (faceHolder.at(i) != NULL)
		{
			delete faceHolder.at(i);
		}
	}
}

//Handles keyboard events
void keyPressed(unsigned char key, int x, int y)
{
	//Quit
	if (key == 'q')
	{
		releaseMemory();
		exit(0);
	}

	//Rotate Up
	if (key == 'w')
	{
		if(cameraPosition[2] < 60)
		cameraPosition[2]++;
	}
	//Rotate down
	if (key == 's')
	{
		if (cameraPosition[2] > -60)
		cameraPosition[2]--;
	}
	//Switch to wireframe
	if (key == 'e')
	{
		viewingMode = true;
	}
	//Switch to faces
	if (key == 'r')
	{
		viewingMode = false;
	}
}

//Set up the inital states of the scene
void init() 
{
	//Turn on depth
	glEnable(GL_DEPTH_TEST);

	//Set Black background
	glClearColor(0.0, 0.0, 0.0, 1.0);
	
	//Set up persepective
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(viewingProjection[0], viewingProjection[1], viewingProjection[2], viewingProjection[3]);

	//Lighting
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GLU_TRUE);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
	glLightfv(GL_LIGHT0, GL_POSITION, startLightPosition);


	//Set Up Camera
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraPosition[0], cameraPosition[1], cameraPosition[2], cameraPosition[3], cameraPosition[4], cameraPosition[5], cameraPosition[6], cameraPosition[7], cameraPosition[8]);
	
}


int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowPosition(80, 80);
	glutInitWindowSize(800, 600);
	glutCreateWindow("Matthew Chiborak CS3388 Assignment 2");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyPressed);
	glutIdleFunc(spinLight);
	
	//Print control instructions
	std::cout << "Press w and s to change camera angle\nPress e and r to switch between wireframe and rendered\nPress q to quit\n";

	init();

	//If cannot read the file, do no display anything
	if (readVertices("vase.txt"))
	{
		generateRotationProfiles();
		triangulatePoints();
		glutMainLoop();
	}

	//Release the dynamic memory used to create the mesh data
	releaseMemory();

	return 0;
}