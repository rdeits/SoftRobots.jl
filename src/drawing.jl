function convert{N, T, PointType}(::Type{Array{T, 2}}, points::Vector{Point{N, PointType}})
	A = Array{T}(N, length(points))
	for i = 1:N
		for j = 1:length(points)
			A[i,j] = points[j][i]
		end
	end
	A
end

function convert{N, T, FaceType, Offset}(::Type{Array{T, 2}}, faces::Vector{GeometryTypes.Face{N, FaceType, Offset}})
    A = Array{T}(N, length(faces))
	for i = 1:N
        for j = 1:length(faces)
            A[i,j] = faces[j][i]
		end
	end
	A
end


function geom_data(mesh::AbstractMesh; color=[1;0;0;0.5])
    msg = lcmdrake[:lcmt_viewer_geometry_data]()
    msg[:type] = msg[:MESH]
    msg[:position] = [0;0;0]
    msg[:quaternion] = [1;0;0;0]
    msg[:color] = color
    msg[:string_data] = ""
    msg[:float_data] = Float64[length(vertices(mesh));
        length(faces(mesh));
        convert(Array{Float64,2}, vertices(mesh))[:];
        convert(Array{Float64,2}, map(f -> convert(Face{3, Int, -1}, f), faces(mesh)))[:];
    ]
    msg[:num_float_data] = length(msg[:float_data])
    msg
end

function link_data(mesh::AbstractMesh, name::AbstractString, robot_num::Integer=1; color=[1;0;0;0.5])
    msg = lcmdrake[:lcmt_viewer_link_data]()
    msg[:name] = name
    msg[:robot_num] = robot_num
    msg[:num_geom] = 1
    push!(msg["geom"], geom_data(mesh, color=color))
    msg
end

function link_data_list(robot::FixedObject, state::FixedObjectState, robot_num::Integer=1)
    lb = Vec(-1., -1, -1)
    ub = Vec(1.5, 1.5, 1.5)
    widths = ub - lb
    bounds = HyperRectangle(lb, widths)
    mesh = barrier_mesh(bounds, state.barrier)
    [link_data(mesh, "barrier", robot_num, color=[0;1;0;0.2])]
end

function link_data_list(robot::SoftRobot, state::SoftRobotState, robot_num::Integer=1)
    [link_data(HomogenousMesh(state.positions, robot.faces), "body", robot_num, color=[1;1;0;0.5])
     link_data(barrier_mesh(robot, state), "barrier", robot_num, color=[1;0;0;0.2])]
 end

function draw(lcm::PyLCM.LCM, world, world_state)
    msg = lcmdrake[:lcmt_viewer_load_robot]()
    for (i, robot) in enumerate(world.objects)
        append!(msg["link"], link_data_list(robot, world_state[robot], i))
    end
    msg[:num_links] = length(msg[:link])
    publish(lcm, "DRAKE_VIEWER_LOAD_ROBOT", msg)
end

function draw(world, world_state)
    lcm = LCM()
    draw(lcm, world, world_state)
    finalize(lcm)
end
