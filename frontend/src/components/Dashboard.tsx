import React, { useEffect, useState, useRef } from "react";
import { useNavigate } from "react-router-dom";
import { FaListAlt, FaArchive, FaUndo } from 'react-icons/fa';
import { useEnvironment } from "../context/EnvironmentContext";
import { useModel } from "../context/ModelContext";
import LoadingSpinner from "./shared/LoadingSpinner";
import { io, Socket } from 'socket.io-client';
import { useSocket } from '../context/SocketContext';
import ChatInputBar from "./shared/ChatInputBar";
import { getWorkflowId, getWorkflowLabelById } from './shared/workflows';
import { FaHammer } from 'react-icons/fa';
import { RiDeleteBin6Line } from 'react-icons/ri';

interface Task {
  task_id: string;
  title: string;
  archived: boolean;
  environment_id: string;
  created_at?: string;
  session_id?: string;
  status?: string; // e.g., 'processing', 'done'
  finished_at?: string;
  workflow_id: number;
}

// const TYPEWRITER_TEXT = "Is it toxic, or just misunderstood? Let's break it down!";
const TYPEWRITER_TEXT = "ToxIndex";
const ENABLE_TYPEWRITER = true; // Set to true to enable typewriter effect
const RESTART_TYPEWRITER = false; // Set to true to restart typewriter after delay
const RESTART_DELAY_MS = 5000; // Delay in ms before restarting typewriter

const Dashboard = () => {
  const [chatInput, setChatInput] = useState("");
  const [error, setError] = useState<string | null>(null);
  const [tasksLoading, setTasksLoading] = useState(false);
  const [activeTasks, setActiveTasks] = useState<Task[]>([]);
  const [archivedTasks, setArchivedTasks] = useState<Task[]>([]);
  const [activeTab, setActiveTab] = useState<'tasks' | 'archive'>('tasks');
  const [uploading, setUploading] = useState(false);
  const { selectedEnv, refetchEnvironments } = useEnvironment();
  const { selectedModel } = useModel();
  const socketRef = useRef<Socket | null>(null);
  const [typedHeading, setTypedHeading] = useState(ENABLE_TYPEWRITER ? "" : TYPEWRITER_TEXT);
  const typewriterTimeoutRef = useRef<number | null>(null);
  const typewriterIntervalRef = useRef<number | null>(null);
  const navigate = useNavigate();
  const prevTaskIdsRef = useRef<string[]>([]);
  const [fileId, setFileId] = useState<string | undefined>(undefined);
  const [fileName, setFileName] = useState<string | undefined>(undefined);

  console.log("Dashboard mounted");

  // Fetch environments when component mounts
  useEffect(() => {
    refetchEnvironments();
  }, [refetchEnvironments]);

  useEffect(() => {
    console.log("selectedEnv", selectedEnv);
    setTasksLoading(true);
    let url = "/api/tasks";
    if (selectedEnv && selectedEnv !== "__add__" && selectedEnv !== "__manage__") {
      url += `?environment_id=${selectedEnv}`;
    }
    fetch(url, { credentials: "include", cache: "no-store" })
      .then(res => res.json())
      .then(data => {
        setActiveTasks(data.active_tasks || []);
        setArchivedTasks(data.archived_tasks || []);
      })
      .catch(err => {
        setActiveTasks([]);
        setArchivedTasks([]);
        setError("Failed to load tasks.");
        console.error("Error loading tasks:", err);
      })
      .finally(() => setTasksLoading(false));
  }, [selectedEnv]);

  useEffect(() => {
    if (!ENABLE_TYPEWRITER) {
      setTypedHeading(TYPEWRITER_TEXT);
      return;
    }
    function startTypewriter() {
      let i = 0;
      setTypedHeading("");
      if (typewriterIntervalRef.current) clearInterval(typewriterIntervalRef.current);
      if (typewriterTimeoutRef.current) clearTimeout(typewriterTimeoutRef.current);
      typewriterIntervalRef.current = setInterval(() => {
        setTypedHeading(TYPEWRITER_TEXT.slice(0, i + 1));
        i++;
        if (i === TYPEWRITER_TEXT.length) {
          if (typewriterIntervalRef.current) clearInterval(typewriterIntervalRef.current);
          if (RESTART_TYPEWRITER) {
            typewriterTimeoutRef.current = setTimeout(() => {
              startTypewriter();
            }, RESTART_DELAY_MS);
          }
        }
      }, 80);
    }
    startTypewriter();
    return () => {
      if (typewriterIntervalRef.current) clearInterval(typewriterIntervalRef.current);
      if (typewriterTimeoutRef.current) clearTimeout(typewriterTimeoutRef.current);
    };
  }, []);

  // Socket.IO setup for real-time task status updates
  useEffect(() => {
    const { socket, isConnected } = useSocket();
    
    if (!socket) return;

    // Only manage rooms if socket is connected
    if (isConnected) {
      // Leave old rooms
      const prevTaskIds = prevTaskIdsRef.current;
      const currentTaskIds = activeTasks.map(t => t.task_id);

      prevTaskIds.forEach(id => {
        if (!currentTaskIds.includes(id)) {
          console.log('[SocketIO] Leaving task room', id);
          socket.emit('leave_task_room', { task_id: id });
        }
      });

      // Join new rooms
      currentTaskIds.forEach(id => {
        if (!prevTaskIds.includes(id)) {
          console.log('[SocketIO] Joining task room', id);
          socket.emit('join_task_room', { task_id: id });
        }
      });

      prevTaskIdsRef.current = currentTaskIds;
    }

    // Global event logger for debugging (only in development)
    if (import.meta.env.DEV) {
      socket.onAny((event: string, ...args: any[]) => {
        console.log('[SocketIO] Event:', event, args);
      });
    }

    // Listen for status updates
    const handler = (data: any) => {
      console.log('[SocketIO] Received task_status_update', data);
      setActiveTasks(prev => {
        const idx = prev.findIndex(t => t.task_id === data.task_id);
        if (idx !== -1) {
          // Update the task in place
          const updated = [...prev];
          updated[idx] = { ...updated[idx], ...data };
          return updated;
        }
        return prev;
      });
    };
    socket.on('task_status_update', handler);

    return () => {
      socket.off('task_status_update', handler);
      if (import.meta.env.DEV) {
        socket.offAny(); // Remove global event logger
      }
      // Don't leave rooms on cleanup - let the server handle it
    };
  }, [activeTasks]);


  const archiveTask = (task_id: string) => {
    fetch(`/api/tasks/${task_id}/archive`, { method: "POST", credentials: "include", cache: "no-store" })
      .then(res => res.json())
      .then(() => {
        setActiveTasks(tasks => tasks.filter(t => t.task_id !== task_id));
        const archived = activeTasks.find(t => t.task_id === task_id);
        if (archived) setArchivedTasks(tasks => [archived, ...tasks]);
      });
  };

  const unarchiveTask = (task_id: string) => {
    fetch(`/api/tasks/${task_id}/unarchive`, { method: "POST", credentials: "include", cache: "no-store" })
      .then(res => res.json())
      .then(() => {
        setArchivedTasks(tasks => tasks.filter(t => t.task_id !== task_id));
        const unarchived = archivedTasks.find(t => t.task_id === task_id);
        if (unarchived) setActiveTasks(tasks => [unarchived, ...tasks]);
      });
  };

  const handleFormSubmit = async (e: React.FormEvent<HTMLFormElement>) => {
    e.preventDefault();
    setUploading(true);
    setError(null);
    if (!chatInput.trim()) {
      setError("Please enter a question.");
      setUploading(false);
      return;
    }
    if (!selectedEnv || selectedEnv === "__add__" || selectedEnv === "__manage__") {
      setError("Please select an environment.");
      setUploading(false);
      return;
    }

    const workflow_id = getWorkflowId(selectedModel);
    if (workflow_id === 4 && !fileId) {
      setError("You must select a file for this workflow.");
      setUploading(false);
      return;
    }
    if (workflow_id === 0) {
      setError("This workflow is not yet supported.");
      setUploading(false);
      return;
    }
    setChatInput(""); // clear input
    try {
      const res = await fetch("/api/tasks", {
        method: "POST",
        credentials: "include",
        cache: "no-store",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          message: chatInput,
          workflow: workflow_id,
          environment_id: selectedEnv,
          file_id: fileId,
        }),
      });
      if (!res.ok) throw new Error("Failed to create task");

      setTasksLoading(true);
      let url = "/api/tasks";
      if (selectedEnv && selectedEnv !== "__add__" && selectedEnv !== "__manage__") {
        url += `?environment_id=${selectedEnv}`;
      }
      const tasksData = await fetch(url, { credentials: "include", cache: "no-store" }).then(r => r.json());
      setActiveTasks(tasksData.active_tasks || []);
      setArchivedTasks(tasksData.archived_tasks || []);
      setTasksLoading(false);
      setFileId(undefined);
      setFileName(undefined);
    } catch (err) {
      setError("Failed to post task.");
    } finally {
      setUploading(false);
    }
  };


  return (
    <div className="min-h-screen flex flex-col flex-1" style={{ background: 'linear-gradient(135deg, #101614 0%, #1a2a1a 60%, #1a2320 100%)', position: 'relative' }}>
      <div className="w-full flex justify-center pt-24 pb-4">
        <h2 className="text-3xl font-normal text-white text-center" style={{ fontFamily: 'inherit', minHeight: 60, letterSpacing: '-0.01em' }}>
          {typedHeading}
          {typedHeading.length < TYPEWRITER_TEXT.length && <span style={{ color: '#6ee7b7', fontWeight: 'bold', marginLeft: 2, animation: 'blink 1s steps(1) infinite' }}>|</span>}
        </h2>
        <style>{`
          @keyframes blink {
            0%, 100% { opacity: 1; }
            50% { opacity: 0; }
          }
        `}</style>
      </div>
      <div className="w-full max-w-2xl mx-auto mt-8">
        <div className="flex space-x-2 mb-4 border-b border-gray-700">
          <button
            className={`px-4 py-2 font-semibold flex items-center gap-2 focus:outline-none transition border-b-4 ${activeTab === 'tasks' ? 'border-white text-white' : 'border-transparent text-gray-400 hover:text-white'}`}
            style={{ background: 'none', borderRadius: 0 }}
            onClick={() => setActiveTab('tasks')}
          >
            <FaListAlt className="inline-block" /> Tasks
          </button>
          <button
            className={`px-4 py-0 font-semibold flex items-center gap-2 focus:outline-none transition border-b-4 ${activeTab === 'archive' ? 'border-white text-white' : 'border-transparent text-gray-400 hover:text-white'}`}
            style={{ background: 'none', borderRadius: 0 }}
            onClick={() => setActiveTab('archive')}
          >
            <FaArchive className="inline-block" /> Archive
          </button>
        </div>
        {activeTab === 'tasks' && (
          tasksLoading ? (
            <div className="text-white">Loading tasks...</div>
          ) : activeTasks.length === 0 ? (
            <div className="text-gray-400">No active tasks.</div>
          ) : (
            <ul className="divide-y">
              {activeTasks.map(task => (
                <li
                  key={task.task_id}
                  className="flex items-center justify-between px-4 py-0 transition-colors duration-150 border-b border-gray-800 hover:border-purple-700 rounded-lg text-base"
                  style={{ minHeight: 44 }}
                >
                  <div
                    className="flex items-center min-w-0 gap-3 cursor-pointer"
                    style={{ flex: '1 1 0%', minWidth: 0 }}
                    onClick={() => (task.status === 'done' || task.status === 'error') && navigate(`/task/${task.task_id}`)}
                    tabIndex={0}
                    role="button"
                    onKeyDown={e => {
                      if ((task.status === 'done' || task.status === 'error') && (e.key === 'Enter' || e.key === ' ')) {
                        navigate(`/task/${task.task_id}`);
                      }
                    }}
                  >
                    <span
                      className="text-white text-left font-medium flex items-center gap-2 flex-nowrap min-w-0"
                      style={{
                        background: 'none',
                        fontSize: '.95rem',
                        overflow: 'hidden',
                        textOverflow: 'ellipsis',
                        whiteSpace: 'nowrap',
                        lineHeight: 1.1,
                        boxShadow: 'none',
                        outline: 'none',
                        border: 'none',
                        maxWidth: 260
                      }}
                      title={task.title}
                    >
                      <span className="truncate min-w-0" style={{maxWidth: 140, display: 'inline-block', verticalAlign: 'middle'}}>{task.title}</span>
                      <span className="ml-2 text-xs text-blue-300 font-mono flex items-center gap-1 flex-nowrap" title="Tool used" style={{whiteSpace: 'nowrap'}}>
                        <FaHammer className="inline-block mr-1 text-blue-400 align-middle text-base" />
                        {getWorkflowLabelById(task.workflow_id)}
                      </span>
                    </span>

                    {task.status === 'error' ? (
                      <span className="ml-2 text-sm text-red-400 font-mono" title="Task failed">
                        &#9888; Error
                      </span>
                    ) : (task.status !== 'done') && (
                      <span className="ml-2 flex items-center">
                        <LoadingSpinner 
                          size="small" 
                          text="" 
                          showTimer={true} 
                          startTime={task.created_at ? new Date(task.created_at).getTime() : undefined} 
                          workflowId={task.workflow_id}
                          status={task.status}
                        />
                        {task.status && (
                          <span className="ml-2 text-xs text-green-300 font-mono" style={{ maxWidth: 120, whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis' }}>
                            {task.status}
                          </span>
                        )}
                      </span>
                    )}

                    {/*
                    {task.created_at && task.status === 'done' && (
                      <span className="ml-2 text-xs text-gray-400 whitespace-nowrap">
                        <span className="font-semibold text-gray-500">Created:</span>

                          {new Date(task.created_at).toLocaleString('en-US', {
                            month: 'short',
                            day: 'numeric',
                            hour: 'numeric',
                            minute: '2-digit',
                            hour12: true,
                          })}

                      </span>
                    )}
                    */}
                    
                    {task.status === 'done' && task.created_at && task.finished_at && (
                      (() => {
                        const start = new Date(task.created_at).getTime();
                        const end = new Date(task.finished_at).getTime();
                        const diff = Math.max(0, end - start);
                        const mins = Math.floor(diff / 60000);
                        const secs = Math.floor((diff % 60000) / 1000);
                        return (
                          <span className="text-xs text-gray-400 ml-2 whitespace-nowrap font-mono" title="Task duration">
                            <span className="font-semibold text-gray-500">Duration:</span>
                            {mins}m {secs}s
                          </span>
                        );
                      })()
                    )}
                    {/* the user does not need to know the session id for now. */}
                    {/*
                    {task.session_id && task.status === 'done' && (
                      <span className="ml-2 px-2 py-0.5 bg-gray-800 rounded text-green-400 font-mono text-xs" title={task.session_id} style={{fontWeight: 500}}>
                        chat {task.session_id.slice(0, 4)}
                      </span>
                    )}
                      */}
                  </div>
                  <div className="w-12" />
                  {task.status === 'done' && (
                    <button
                      onClick={e => {
                        e.stopPropagation();
                        archiveTask(task.task_id);
                      }}
                      className="ml-2 flex items-center justify-center w-7 h-7 rounded-full !bg-gray-700 hover:bg-gray-600 active:bg-gray-800 text-white transition border-none shadow-sm focus:outline-none focus:ring-2 focus:ring-gray-400 hover:ring-2 hover:ring-gray-400"
                      style={{ padding: 0, borderRadius: '50%', fontSize: '1.1rem', border: 'none', outline: 'none', boxShadow: 'none', cursor: 'pointer' }}
                      title="Archive task"
                    >
                      <FaArchive />
                    </button>
                  )}
                </li>
              ))}
            </ul>
          )
        )}
        {activeTab === 'archive' && (
          tasksLoading ? (
            <div className="text-white">Loading tasks...</div>
          ) : archivedTasks.length === 0 ? (
            <div className="text-gray-400">No archived tasks.</div>
          ) : (
            <ul className="divide-y">
              {archivedTasks.map(task => (
                <li
                  key={task.task_id}
                  className="flex items-center justify-between px-4 py-3 border-b border-gray-800 hover:border-purple-800 rounded"
                  style={{ minHeight: 44 }}
                >
                  <div
                    className="flex items-center min-w-0 gap-3 cursor-pointer"
                    style={{ flex: '1 1 0%', minWidth: 0 }}
                    onClick={() => navigate(`/task/${task.task_id}`)}
                    tabIndex={0}
                    role="button"
                    onKeyDown={e => {
                      if (e.key === 'Enter' || e.key === ' ') navigate(`/task/${task.task_id}`);
                    }}
                  >
                    <span className="text-gray-300 flex items-center gap-2 flex-nowrap min-w-0">
                      <span className="truncate min-w-0" style={{maxWidth: 140, display: 'inline-block', verticalAlign: 'middle'}}>{task.title}</span>
                      <span className="ml-2 text-xs text-blue-300 font-mono flex items-center gap-1 flex-nowrap" title="Tool used" style={{whiteSpace: 'nowrap'}}>
                        <FaHammer className="inline-block mr-1 text-blue-400 align-middle text-base" />
                        {getWorkflowLabelById(task.workflow_id)}
                      </span>
                      <span className="ml-2 text-xs text-gray-400">
                        <span className="font-bold text-gray-500">Archived - </span>{' task initially created at: '}
                        {task.created_at &&
                          new Date(task.created_at).toLocaleString(undefined, {
                            month: 'short', day: 'numeric'
                          })}
                      </span>
                    </span>
                  </div>
                  <div className="w-12" />
                  <button
                    onClick={e => {
                      e.stopPropagation();
                      unarchiveTask(task.task_id);
                    }}
                    className="ml-1 px-2 rounded bg-green-700 text-white hover:bg-green-800 text-xs flex items-center justify-center"
                    style={{ marginTop: 0, marginBottom: 0, height: 'auto', minHeight: 0, paddingTop: 0, paddingBottom: 0, background: 'none', border: 'none', outline: 'none', boxShadow: 'none', cursor: 'pointer' }}
                    title="Unarchive task"
                  >
                    <FaUndo className="mr-1" />
                  </button>
                </li>
              ))}
            </ul>
          )
        )}
      </div>
      <div className="flex-1 flex flex-col justify-end">
        <div className="w-full flex flex-col items-center pb-8">
          <div style={{ width: '100%', maxWidth: 800 }}>
            {fileId && (
              <div className="bg-gray-900/10 text-green-200 px-3 py-0 rounded-full text-sm font-medium flex items-center gap-1 mb-2" style={{ minHeight: 40 }}>
                <span className="font-mono font-semibold text-base z-13">Selected file:</span>
                <span className="font-mono truncate max-w-xs text-base z-13">
                  {fileName
                    ? (() => {
                        if (fileName.length > 24) {
                          const extMatch = fileName.match(/(\.[^./\\]+)$/);
                          const ext = extMatch ? extMatch[1] : '';
                          const base = fileName.slice(0, 20);
                          return base + '...' + ext;
                        } else {
                          return fileName;
                        }
                      })()
                    : fileId.slice(0, 10) + '...'}
                </span>
                <button
                  onClick={() => { setFileId(undefined); setFileName(undefined); }}
                  className="ml-0 text-red-300 hover:text-red-500 p-1 bg-transparent border-none outline-none focus:outline-none"
                  style={{ background: 'none', border: 'none', outline: 'none', boxShadow: 'none', padding: 4, margin: 0 }}
                  title="Clear file"
                >
                  <RiDeleteBin6Line className="w-6 h-6" />
                </button>
              </div>
            )}
            <ChatInputBar
              value={chatInput}
              onChange={setChatInput}
              onSubmit={handleFormSubmit}
              uploading={uploading}
              error={error}
              onFilePick={(id: string, name?: string) => {
                setFileId(id);
                setFileName(name);
              }}
            />
          </div>
        </div>
      </div>
    </div>
  );
};

export default Dashboard;