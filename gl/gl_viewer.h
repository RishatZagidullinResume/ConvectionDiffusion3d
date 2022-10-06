#pragma once
#include <GLFW/glfw3.h>
#include <glm/gtc/type_ptr.hpp>
#include "shader.h"

namespace custom_gl
{
    class gl_viewer
    {
    private:
        short SCR_WIDTH;
        short SCR_HEIGHT;
        Shader *fill_shader;
        Shader *line_shader;
    public:
        GLFWwindow *window;
        unsigned int VAO, VBO, EBO, VBO2, EBO2;
        void view(size_t const &size, size_t const & shape_size, int time);
        void preset();
        void buffer_ebo(size_t const &size,
                        unsigned int const * ebo = nullptr);
        void buffer_shape_ebo(size_t const &size,
                        unsigned int const * ebo = nullptr);
        void buffer_vbo(size_t const &size,
                        float const * vbo = nullptr);
        void buffer_shape_vbo(size_t const &size,
                        float const * vbo = nullptr);
        gl_viewer(short SCR_WIDTH, short SCR_HEIGHT, bool hide_window);
        ~gl_viewer();
    };

    gl_viewer::gl_viewer(short SCR_WIDTH, short SCR_HEIGHT, bool hide_window=false)
    {
        this->SCR_WIDTH = SCR_WIDTH;
        this->SCR_HEIGHT = SCR_HEIGHT;
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        if (hide_window)
            glfwWindowHint(GLFW_VISIBLE, GL_FALSE);

        window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT,
                                  "", NULL, NULL);
        glfwMakeContextCurrent(window);
        gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
        const GLubyte* vendor = glGetString(GL_VENDOR);
        const GLubyte* renderer = glGetString(GL_RENDERER);
        std::cout << vendor << "  " << renderer << "\n";

        fill_shader = new Shader("shaders/fill_shader.vs",
                                 "shaders/fill_shader.fs");
        line_shader = new Shader("shaders/line_shader.vs",
                                 "shaders/line_shader.fs");
        glEnable(GL_DEPTH_TEST);
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);
        glGenBuffers(1, &VBO2);
        glGenBuffers(1, &EBO2);
        glBindVertexArray(VAO);
    }

    void gl_viewer::view(size_t const &size, size_t const &shape_size, int time = -1)
    {
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, true);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        float aaa[16] = {
            1.810660, 0.000000, 0.000000, 0.000000,
            0.000000, 2.414213, 0.000000, 0.000000,
            0.000000, 0.000000, -1.002002, -1.000000,
            0.000000, 0.000000, -0.200200, 0.000000
        };
        glm::mat4 projection;
        memcpy( glm::value_ptr(projection), aaa, sizeof(aaa));
        float bbb[16] = {
            0.754713, -0.251062, 0.606116, 0.000000,
            -0.656055, -0.288817, 0.697264, 0.000000,
            0.000000, -0.923879, -0.382684, 0.000000,
            0.038596, -0.097829, -5.060080, 1.000000
        };
        glm::mat4 view;
        memcpy( glm::value_ptr(view), bbb, sizeof(bbb));
        
        fill_shader->use();
        fill_shader->setMat4("projection", projection);
        fill_shader->setMat4("view", view);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0, 1.0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        preset(); 

        glDrawElements(GL_TRIANGLES, size, GL_UNSIGNED_INT, 0);

        glBindBuffer(GL_ARRAY_BUFFER, VBO2);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO2);
        preset(); 

        glDrawElements(GL_TRIANGLES, shape_size, GL_UNSIGNED_INT, 0);

        glDisable(GL_POLYGON_OFFSET_FILL);
        line_shader->use();
        line_shader->setMat4("projection", projection);
        line_shader->setMat4("view", view);
        glLineWidth(5.0);
        glPointSize(5.0);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        preset(); 

        glDrawElements(GL_LINES, size, GL_UNSIGNED_INT, 0);

        glBindBuffer(GL_ARRAY_BUFFER, VBO2);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO2);
        preset(); 

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawElements(GL_TRIANGLES, shape_size, GL_UNSIGNED_INT, 0);


        glfwSwapBuffers(window);
        glfwPollEvents();

        if (time != -1)
        {
            int * buffer = new int [SCR_WIDTH * SCR_HEIGHT * 3];
            glReadPixels(0, 0, SCR_WIDTH, SCR_HEIGHT, GL_BGR,
                         GL_UNSIGNED_BYTE, buffer);
            std::string filename = std::string("imgs/data_")
                    + std::to_string(time) + std::string(".tga");
            FILE *out = std::fopen(filename.c_str(), "w");
            short head[] {0,2,0,0,0,0,SCR_WIDTH,SCR_HEIGHT,24};
            fwrite(&head, sizeof(head), 1, out);
            fwrite(buffer, 3*SCR_WIDTH*SCR_HEIGHT, 1, out);
            fclose(out);
            delete [] buffer;
        }
        return;
    }

    void gl_viewer::preset()
    {
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
                              4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer( 1, 1, GL_FLOAT, GL_FALSE,
                              4 * sizeof(float),
                              (void*)(3*sizeof(float)) );
        glEnableVertexAttribArray(1);
        return;
    }

    void gl_viewer::buffer_ebo(size_t const &size,
                               unsigned int const * ebo)
    {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                     sizeof(unsigned int)*size, ebo, GL_STATIC_DRAW);
        return;
    }
    
    void gl_viewer::buffer_shape_ebo(size_t const &size,
                               unsigned int const * ebo)
    {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO2);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                     sizeof(unsigned int)*size, ebo, GL_STATIC_DRAW);
        return;
    }

    void gl_viewer::buffer_vbo(size_t const &size, float const * vbo)
    {
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float)*size,
                     vbo, GL_DYNAMIC_DRAW);
        return;
    }
    
    void gl_viewer::buffer_shape_vbo(size_t const &size, float const * vbo)
    {
        glBindBuffer(GL_ARRAY_BUFFER, VBO2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float)*size,
                     vbo, GL_DYNAMIC_DRAW);
        return;
    }

    gl_viewer::~gl_viewer()
    {
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &VBO);
        glDeleteBuffers(1, &EBO);
        glfwTerminate();
    }
}





