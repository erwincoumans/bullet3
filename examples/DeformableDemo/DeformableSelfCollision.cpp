/*
 Bullet Continuous Collision Detection and Physics Library
 Copyright (c) 2019 Google Inc. http://bulletphysics.org
 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from the use of this software.
 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it freely,
 subject to the following restrictions:
 1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
 2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
 3. This notice may not be removed or altered from any source distribution.
 */
#include "DeformableSelfCollision.h"
///btBulletDynamicsCommon.h is the main Bullet include file, contains most common include files.
#include "btBulletDynamicsCommon.h"
#include "BulletSoftBody/btDeformableMultiBodyDynamicsWorld.h"
#include "BulletSoftBody/btSoftBody.h"
#include "BulletSoftBody/btSoftBodyHelpers.h"
#include "BulletSoftBody/btDeformableBodySolver.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraintSolver.h"
#include "../CommonInterfaces/CommonParameterInterface.h"
#include <stdio.h>  //printf debugging

#include "../CommonInterfaces/CommonDeformableBodyBase.h"
#include "../Utils/b3ResourcePath.h"
///////////////////////////////////////
#include <OpenTissue/core/containers/t4mesh/util/t4mesh_block_generator.h>
#include <OpenTissue/core/math/math_basic_types.h>
#include <OpenTissue/dynamics/fem/fem.h>
#include "MyFemMesh.h"
typedef OpenTissue::math::BasicMathTypes<double, size_t>    math_types;
typedef OpenTissue::fem::Mesh<math_types>                  mesh_type;
typedef math_types::vector3_type                           vector3_type;
typedef math_types::real_type                              real_type;
bool m_fixNodes = true;
static void toggleFixedNodes(int buttonId, bool buttonState, void* userPointer)
{
    m_fixNodes = !m_fixNodes;
    // todo(thomas) is there a way to get a toggle button with changing text?
    b3Printf("switched fixed nodes %s", m_fixNodes ? "on" : "off");
}


struct world_point_container
{
    typedef vector3_type value_type;

    world_point_container(mesh_type* mesh)
        : m_mesh(mesh)
    {}

    vector3_type& operator[] (unsigned int const& idx)
    {
        return m_mesh->m_nodes[idx].m_coord;
    }

    vector3_type const& operator[] (unsigned int const& idx)const
    {
        return m_mesh->m_nodes[idx].m_coord;
    }

    void clear() {}
    unsigned int size()const { return m_mesh->m_nodes.size(); }
    void resize(unsigned int) {}

    mesh_type* m_mesh;
};
///////////////////////////////////////

///The DeformableSelfCollision shows deformable self collisions
class DeformableSelfCollision : public CommonDeformableBodyBase
{
 ///////////////////////////////////////
    real_type		m_young;// = 500000;
    real_type		m_poisson;// = 0.33;
    real_type		m_density;// = 1000;

    //--- infinite m_c_yield plasticity settings means that plasticity is turned off
    real_type		m_c_yield;// = .04;  //--- should be less than maximum expected elastic strain in order to see effect (works as a minimum).
    real_type		m_c_creep;// = .20;  //--- controls how fast the plasticity effect occurs (it is a rate-like control).
    real_type		m_c_max;// = 0.2;    //--- This is maximum allowed plasticity strain (works as a maximum).

    bool m_stiffness_warp_on;
    bool m_collideGroundPlane;
    
    real_type m_gravity;
    mesh_type		m_mesh1;
///////////////////////////////////////
public:
    DeformableSelfCollision(struct GUIHelperInterface* helper)
    : CommonDeformableBodyBase(helper),
        m_stiffness_warp_on(true),
        m_collideGroundPlane(true),
        m_gravity(29.81)
    {
    }
    virtual ~DeformableSelfCollision()
    {
    }
    void initPhysics();
    
    void exitPhysics();
    
    void resetCamera()
    {
        float dist = 3.2;
        float pitch = -17;
        float yaw = -50;
        float targetPos[3] = {0.11, 0.60, 0.16};
        m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
    }
    
    void stepSimulation(float deltaTime)
    {
        float internalTimeStep = 1. / 240.f;
        //m_dynamicsWorld->stepSimulation(deltaTime, 4, internalTimeStep);


///////////////////////////////////////



        double damp = 2.f;
        double dt = 0.001;

        // Calculate external forces acting on the model.
        // The external forces are stored in each vertex and consists of a downward gravity.
        // If a vertex is being dragged by the user, its velocity vector is added as an external force.

        for (int n = 0; n < m_mesh1.m_nodes.size(); n++)
        {

            if (m_fixNodes)
            {
                if (m_mesh1.m_nodes[n].m_model_coord(0) < 0.01)
                {
                    m_mesh1.m_nodes[n].m_fixed = true;
                }
            }
            else
            {
                if (m_mesh1.m_nodes[n].m_model_coord(0) < 0.01)
                {
                    m_mesh1.m_nodes[n].m_fixed = false;
                }
            }

            if (m_collideGroundPlane && m_mesh1.m_nodes[n].m_coord(1) < 0.f)
            {
                float depth = -m_mesh1.m_nodes[n].m_coord(1);
                if (depth > 0.1)
                    depth = 0.1;

                //#define USE_FORCE
                //#ifdef USE_FORCE
                //			m_mesh1.m_nodes[n].m_f_external = vector3_type(0.0, depth*100000, 0.0);
                //#else
                            //m_mesh1.m_nodes[n].m_coord(1) += 0.95*depth;

                m_mesh1.m_nodes[n].m_f_external = vector3_type(0.0,  depth * 10000, 0.0 );
                //			float dt = 1./60.f;

                if (m_mesh1.m_nodes[n].m_velocity(1) < 0.f)
                {
                    m_mesh1.m_nodes[n].m_velocity(1) = 0.f;
                }

                m_mesh1.m_nodes[n].m_velocity(0) = 0.f;
                m_mesh1.m_nodes[n].m_velocity(2) = 0.f;

                //	m_mesh1.m_nodes[n].m_f_external(1) -= 0.01f*m_mesh1.m_nodes[n].m_velocity(1)/dt;///dt;// depth*100;



      //#endif
            }
            else
            {
                //normal gravity
                m_mesh1.m_nodes[n].m_f_external = vector3_type(0.0, -(m_mesh1.m_nodes[n].m_mass * m_gravity), 0.0);
            }
        }


        /*for (int n=0;n<m_mesh2.m_nodes.size();n++)
        {
            m_mesh2.m_nodes[n].m_f_external = vector3_type(0.0, -(m_mesh2.m_nodes[n].m_mass * m_gravity), 0.0);
        }
        */


        //OpenTissue::fem::simulate(m_mesh1,0.01,m_stiffness_warp_on,1.0,20,20);
        {
            B3_PROFILE("fem::simulate");
            OpenTissue::fem::simulate(m_mesh1, dt, m_stiffness_warp_on, damp);//,0.1,20,20);//,1.0,20,20);
        }

        struct tmp_point
        {
            float m_x, m_y, m_z;
            tmp_point(float x, float y, float z):
                m_x(x), m_y(y), m_z(z)

            {
            }
            
        };
        std::vector<tmp_point> positions;
        std::vector<unsigned int> indices;
        {
            B3_PROFILE("fem::render");
            world_point_container point_wrapper(&m_mesh1);
            for (int tet = 0; tet < point_wrapper.m_mesh->m_tetrahedra.size(); tet++)
            {
                auto tetra = point_wrapper.m_mesh->m_tetrahedra[tet];
                {
                    //T.scale(scale);
                    //DrawTetrahedron(T, wireframe);
                    for (int i = 0; i < 4; i++)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            auto p0 = point_wrapper[tetra.m_nodes[i]];
                            auto p1 = point_wrapper[tetra.m_nodes[j]];
                            tmp_point from(p0.x, p0.y, p0.z);
                            tmp_point to(p1.x, p1.y, p1.z);
                            indices.push_back(positions.size());
                            positions.push_back(from);
                            indices.push_back(positions.size());
                            positions.push_back(to);
                        }
                    }

                }
            }
        }

        {
            float color[4] = { 0,0,0,1 };
            float lineWidthIn = 2;
            m_guiHelper->getRenderInterface()->drawLines(&positions[0].m_x,
                color, positions.size(),
                sizeof(tmp_point),
                &indices[0], indices.size(), lineWidthIn);
            /*app.m_renderer->draw_lines(&positions[0],
                color,
                positions.size(), sizeof(btVector3),
                &indices[0],
                indices.size(), lineWidthIn);
                */
        }

///////////////////////////////////////


    }
    
    void addCloth(btVector3 origin);
    
    virtual void renderScene()
    {
        CommonDeformableBodyBase::renderScene();
    }
};

void DeformableSelfCollision::initPhysics()
{
    m_guiHelper->setUpAxis(1);
    ///collision configuration contains default setup for memory, collision setup
    m_collisionConfiguration = new btSoftBodyRigidBodyCollisionConfiguration();
    
    ///use the default collision dispatcher. For parallel processing you can use a diffent dispatcher (see Extras/BulletMultiThreaded)
    m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);
    
    m_broadphase = new btDbvtBroadphase();
    btDeformableBodySolver* deformableBodySolver = new btDeformableBodySolver();
    
    ///the default constraint solver. For parallel processing you can use a different solver (see Extras/BulletMultiThreaded)
    btDeformableMultiBodyConstraintSolver* sol = new btDeformableMultiBodyConstraintSolver();
    sol->setDeformableSolver(deformableBodySolver);
    m_solver = sol;
    
    m_dynamicsWorld = new btDeformableMultiBodyDynamicsWorld(m_dispatcher, m_broadphase, sol, m_collisionConfiguration, deformableBodySolver);
    //    deformableBodySolver->setWorld(getDeformableDynamicsWorld());
    //    m_dynamicsWorld->getSolverInfo().m_singleAxisDeformableThreshold = 0.f;//faster but lower quality
    btVector3 gravity = btVector3(0, -9.8, 0);
    m_dynamicsWorld->setGravity(gravity);
    getDeformableDynamicsWorld()->getWorldInfo().m_gravity = gravity;
    getDeformableDynamicsWorld()->getWorldInfo().m_sparsesdf.setDefaultVoxelsz(0.25);
    
    //    getDeformableDynamicsWorld()->before_solver_callbacks.push_back(dynamics);
    m_guiHelper->createPhysicsDebugDrawer(m_dynamicsWorld);
    
    {
        ///create a ground
        btCollisionShape* groundShape = new btBoxShape(btVector3(btScalar(150.), btScalar(2.5), btScalar(150.)));
        groundShape->setMargin(0.0001);
        m_collisionShapes.push_back(groundShape);
        
        btTransform groundTransform;
        groundTransform.setIdentity();
        groundTransform.setOrigin(btVector3(0, -3.5, 0));
        groundTransform.setRotation(btQuaternion(btVector3(1, 0, 0), SIMD_PI * 0));
        //We can also use DemoApplication::localCreateRigidBody, but for clarity it is provided here:
        btScalar mass(0.);
        
        //rigidbody is dynamic if and only if mass is non zero, otherwise static
        bool isDynamic = (mass != 0.f);
        
        btVector3 localInertia(0, 0, 0);
        if (isDynamic)
            groundShape->calculateLocalInertia(mass, localInertia);
        
        //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
        btDefaultMotionState* myMotionState = new btDefaultMotionState(groundTransform);
        btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, groundShape, localInertia);
        btRigidBody* body = new btRigidBody(rbInfo);
        body->setFriction(4);
        
        //add the ground to the dynamics world
        m_dynamicsWorld->addRigidBody(body);
    }
    addCloth(btVector3(0, -0.2, 0));
    addCloth(btVector3(0, -0.1, 0));
    getDeformableDynamicsWorld()->setImplicit(false);
    getDeformableDynamicsWorld()->setLineSearch(false);
    m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);


///////////////////////////////////////
    
    ButtonParams button("toggle fixed nodes", 0, true);
    button.m_callback = toggleFixedNodes;
    this->m_guiHelper->getParameterInterface()->registerButtonParameter(button);

    world_point_container point_wrapper(&m_mesh1);
    OpenTissue::t4mesh::generate_blocks(10, 3, 3, 0.1, 0.1, 0.1, m_mesh1);

    for (int n = 0; n < m_mesh1.m_nodes.size(); n++)
    {
        m_mesh1.m_nodes[n].m_coord(1) += 2.0f;

        m_mesh1.m_nodes[n].m_model_coord = m_mesh1.m_nodes[n].m_coord;

    }

    


    m_young = 1000000;//47863;//100000;
    m_poisson = 0.33;//33;
  //m_poisson = 0.48;
    m_density = 1054.00;//1000;

    //--- infinite m_c_yield plasticity settings means that plasticity is turned off
    m_c_yield = 0.03;//.04;  //--- should be less than maximum expected elastic strain in order to see effect (works as a minimum).
    m_c_creep = 0.20;//.20;  //--- controls how fast the plasticity effect occurs (it is a rate-like control).
    m_c_max = 1e30f;//0.2;    //--- This is maximum allowed plasticity strain (works as a maximum).


    OpenTissue::fem::init(m_mesh1, m_young, m_poisson, m_density, m_c_yield, m_c_creep, m_c_max);


///////////////////////////////////////

}
void DeformableSelfCollision::addCloth(btVector3 origin)
// create a piece of cloth
{
    const btScalar s = 0.6;
    const btScalar h = 0;
    
    btSoftBody* psb = btSoftBodyHelpers::CreatePatch(getDeformableDynamicsWorld()->getWorldInfo(), btVector3(-s, h, -2*s),
                                                     btVector3(+s, h, -2*s),
                                                     btVector3(-s, h, +2*s),
                                                     btVector3(+s, h, +2*s),
                                                     15,30,
//                                                     4,4,
                                                     0, true, 0.0);

    
    psb->getCollisionShape()->setMargin(0.02);
    psb->generateBendingConstraints(2);
    psb->setTotalMass(.5);
    psb->m_cfg.kKHR = 1; // collision hardness with kinematic objects
    psb->m_cfg.kCHR = 1; // collision hardness with rigid body
    psb->m_cfg.kDF = 0.1;
//    psb->rotate(btQuaternion(0, SIMD_PI / 2, 0));
    btTransform clothTransform;
    clothTransform.setIdentity();
    clothTransform.setOrigin(btVector3(0,0.2,0)+origin);
    psb->transform(clothTransform);
    psb->m_cfg.collisions = btSoftBody::fCollision::SDF_RD;
    psb->m_cfg.collisions |= btSoftBody::fCollision::SDF_RDF;
    psb->m_cfg.collisions |= btSoftBody::fCollision::VF_DD;
    getDeformableDynamicsWorld()->addSoftBody(psb);
    psb->setSelfCollision(true);
    
    btDeformableMassSpringForce* mass_spring = new btDeformableMassSpringForce(2,0.2, true);
    psb->setSpringStiffness(4);
    getDeformableDynamicsWorld()->addForce(psb, mass_spring);
    m_forces.push_back(mass_spring);
    btVector3 gravity = btVector3(0, -9.8, 0);
    btDeformableGravityForce* gravity_force =  new btDeformableGravityForce(gravity);
    getDeformableDynamicsWorld()->addForce(psb, gravity_force);
    m_forces.push_back(gravity_force);
}

void DeformableSelfCollision::exitPhysics()
{
    //cleanup in the reverse order of creation/initialization
    removePickingConstraint();
    //remove the rigidbodies from the dynamics world and delete them
    int i;
    for (i = m_dynamicsWorld->getNumCollisionObjects() - 1; i >= 0; i--)
    {
        btCollisionObject* obj = m_dynamicsWorld->getCollisionObjectArray()[i];
        btRigidBody* body = btRigidBody::upcast(obj);
        if (body && body->getMotionState())
        {
            delete body->getMotionState();
        }
        m_dynamicsWorld->removeCollisionObject(obj);
        delete obj;
    }
    // delete forces
    for (int j = 0; j < m_forces.size(); j++)
    {
        btDeformableLagrangianForce* force = m_forces[j];
        delete force;
    }
    m_forces.clear();
    //delete collision shapes
    for (int j = 0; j < m_collisionShapes.size(); j++)
    {
        btCollisionShape* shape = m_collisionShapes[j];
        delete shape;
    }
    m_collisionShapes.clear();
    
    delete m_dynamicsWorld;
    
    delete m_solver;
    
    delete m_broadphase;
    
    delete m_dispatcher;
    
    delete m_collisionConfiguration;
}



class CommonExampleInterface* DeformableSelfCollisionCreateFunc(struct CommonExampleOptions& options)
{
    return new DeformableSelfCollision(options.m_guiHelper);
}


