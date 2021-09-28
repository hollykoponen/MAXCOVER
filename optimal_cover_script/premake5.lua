workspace "optimal_covers"
   startproject "optimal_covers"
   configurations { "Debug", "Release" }
--   cleancommands { "{RMDIR} ".. path.join("include", "RMQ") } -- Doesn't work...???

project "optimal_covers"
   kind "ConsoleApp"
   language "C++"
   architecture "x86_64"

   targetdir ( path.join("bin", "%{cfg.buildcfg}") )
   objdir ( path.join("obj", "%{cfg.buildcfg}") )

   dependson {
      "RMQ",
      "lcpdc",
      "sais"
   }
--   prebuildcommands {
--      "{RMDIR} " .. path.join("include", "lcpdc"),
--      "{COPY} " .. path.join("dependencies", "lcpdc", "include") .. " " .. path.join("include", "lcpdc"),
--      "{RMDIR} " .. path.join("include", "lcpdc", "RMQ"),
--      "{RMDIR} " .. path.join("include", "RMQ"),
--      "{COPY} " .. path.join("dependencies", "RMQ", "include") .. " " .. path.join("include", "RMQ"),
--      "{RMDIR} " .. path.join("include", "sais"),
--      "{COPY} " .. path.join("dependencies", "sais-lite-2.4.1.2", "include") .. " " .. path.join("include", "sais")
--   }

   libdirs { "lib" }
   includedirs { "include" }

   files { path.join("src", "**.cpp") }

--   defines { "__cplusplus" }

   filter "configurations:Debug"
      defines { "DEBUG" }
      symbols "On"

   filter "configurations:Release"
      defines { "NDEBUG" }
      optimize "On"

   filter "system:linux"
      links { "RMQ", "lcpdc", "sais"}
    --   buildoptions "-std=c++2a"

   filter "system:windows"
      defines { "WIN" }
      links { "RMQ", "lcpdc", "sais" }
    --   cppdialect "C++17"


project "RMQ"
   kind "StaticLib"
   language "C++"
   architecture "x86_64"
   buildoptions { "-Wno-narrowing" }

   targetdir "lib"
   objdir ( path.join("obj", "%{cfg.buildcfg}") )

   includedirs { path.join("dependencies", "RMQ", "include") }

   files { path.join("dependencies", "RMQ", "src", "**.cpp") }

   filter "configurations:Debug"
      defines { "DEBUG" }
      symbols "On"

   filter "configurations:Release"
      defines { "NDEBUG" }
      optimize "On"

   filter "system:linux"
      buildoptions { "-Wno-narrowing" }
--      buildoptions "-std=c++2a"

   filter "system:windows"
      defines { "WIN" }
--      cppdialect "C++17"


project "lcpdc"
   kind "StaticLib"
   language "C++"
   architecture "x86_64"

   dependson {
      "RMQ"
   }
--   prebuildcommands {
--      "{RMDIR} " .. path.join("dependencies", "lcpdc", "include", "RMQ"),
--      "{COPY} " .. path.join("dependencies", "RMQ", "include") .. " " .. path.join("dependencies", "lcpdc", "include", "RMQ")
--   }

   targetdir "lib"
   objdir ( path.join("obj", "%{cfg.buildcfg}") )

   includedirs { path.join("dependencies", "lcpdc", "include") }

   files { path.join("dependencies", "lcpdc", "src", "**.cpp") }

   filter "configurations:Debug"
      defines { "DEBUG" }
      symbols "On"

   filter "configurations:Release"
      defines { "NDEBUG" }
      optimize "On"

   filter "system:linux"
      links { "RMQ" }
--      buildoptions "-std=c++2a"

   filter "system:windows"
      links { "RMQ" }
      defines { "WIN" }
--      cppdialect "C++17"



project "sais"
   kind "StaticLib"
   language "C"
   architecture "x86_64"

   targetdir "lib"
   objdir ( path.join("obj", "%{cfg.buildcfg}") )

   includedirs { path.join("dependencies", "sais-lite-2.4.1.2", "include") }

   files { path.join("dependencies", "sais-lite-2.4.1.2", "src", "**.c") }


   filter "configurations:Debug"
      defines { "DEBUG" }
      symbols "On"

   filter "configurations:Release"
      defines { "NDEBUG" }
      optimize "On"

   filter "system:linux"
   --      buildoptions "-std=c++2a"

   filter "system:windows"
      defines { "WIN" }
--      cppdialect "C++17"
