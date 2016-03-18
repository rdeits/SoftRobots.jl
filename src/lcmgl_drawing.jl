
function draw(lcmgl::LCMGLClient, mesh::AbstractMesh; face_color=[.7, .7, .2, .4], draw_normals=false)
    color(lcmgl, face_color...)
    begin_mode(lcmgl, LCMGL.TRIANGLES)
    verts = vertices(mesh)
    for face in faces(mesh)
        for point in verts[face][end:-1:1]
            vertex(lcmgl, point[1], point[2], point[3])
        end
    end
    end_mode(lcmgl)
    color(lcmgl, 0., 0., 0., 1.0)
    begin_mode(lcmgl, LCMGL.LINES)
    for face in faces(mesh)
        for i = 1:length(verts[face])
            vertex(lcmgl, verts[face][i]...)
            j = i + 1
            if i == length(verts[face])
                j = 1
            end
            vertex(lcmgl, verts[face][j]...)
        end
    end
    if draw_normals
      normals = decompose(Normal{3, Float64}, mesh)
      for i = 1:length(verts)
        vertex(lcmgl, verts[i]...)
        vertex(lcmgl, (verts[i] + Point{3, Float64}(normals[i]) * 0.05)...)
      end
    end
    end_mode(lcmgl)
end

function draw(mesh::AbstractMesh)
    LCMGLClient("soft_robot") do lcmgl
      draw(lcmgl, mesh)
      switch_buffer(lcmgl)
    end
end

function draw(robot::SoftRobot, state::SoftRobotState)
  LCMGLClient("soft_robot") do lcmgl
    draw(lcmgl, robot, state)
    switch_buffer(lcmgl)
  end
end



function draw(lcmgl::LCMGLClient, robot::SoftRobot, state::SoftRobotState)
  body_mesh = HomogenousMesh(state.positions, robot.faces)
  draw(lcmgl, body_mesh)

  LCMGLClient("$(lcmgl.name)_barrier") do lcmgl
    draw(lcmgl, barrier_mesh(robot, state), face_color=[.8, .1, .1, .3])
    switch_buffer(lcmgl)
  end
end
