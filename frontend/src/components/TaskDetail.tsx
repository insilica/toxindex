import React, { useEffect, useState, useCallback } from "react";
import { useParams, useNavigate } from "react-router-dom";
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import HomeButton from "./shared/HomeButton";
import FilePreviewModal from './FilePreviewModal';

interface Task {
  task_id: string;
  title: string;
  archived?: boolean;
  environment_id?: string;
  created_at?: string;
  finished_at?: string;
  session_id?: string;
  user_id?: string;
  workflow_id?: number;
}

const getDuration = (start?: string, end?: string) => {
  if (!start || !end) return null;
  const startTime = new Date(start).getTime();
  const endTime = new Date(end).getTime();
  const diff = Math.max(0, endTime - startTime);
  const mins = Math.floor(diff / 60000);
  const secs = Math.floor((diff % 60000) / 1000);
  return `${mins}m ${secs}s`;
};

const TaskDetail: React.FC = () => {
  const { task_id } = useParams<{ task_id: string }>();
  const [task, setTask] = useState<Task | null>(null);
  const [loading, setLoading] = useState(true);
  const [assistantMessage, setAssistantMessage] = useState<string | null>(null);
  const [messageLoading, setMessageLoading] = useState(false);
  const navigate = useNavigate();

  // Tab state
  // const [activeTab, setActiveTab] = useState<'MDrender' | 'MDraw' | 'JsonSchema'>('MDrender');
  const [toxicitySchema, setToxicitySchema] = useState<any>(null);

  const [taskFiles, setTaskFiles] = useState<any[]>([]);
  const [envFiles, setEnvFiles] = useState<any[]>([]);

  const [previewFileId, setPreviewFileId] = useState<string | null>(null);
  const [previewOpen, setPreviewOpen] = useState(false);

  const fetchTask = useCallback(async () => {
    if (!task_id) return;
    setLoading(true);
    try {
      const res = await fetch(`/api/tasks/${task_id}`, { credentials: "include" });
      if (!res.ok) {
        setTask(null);
        setLoading(false);
        return;
      }
      const data = await res.json();
      if (data.error) {
        setTask(null);
      } else {
        setTask(data);
      }
    } catch (e) {
      setTask(null);
    } finally {
      setLoading(false);
    }
  }, [task_id]);

  // Fetch assistant message for this task
  useEffect(() => {
    const fetchAssistantMessage = async () => {
      if (!task_id) return;
      setMessageLoading(true);
      try {
        const res = await fetch(`/api/tasks/${task_id}/messages`, { credentials: "include" });
        if (!res.ok) throw new Error('Failed to fetch messages');
        const data = await res.json();
        if (data.messages && Array.isArray(data.messages)) {
          const assistantMsg = data.messages.find((m: any) => m.role === "assistant");
          if (assistantMsg) {
            setAssistantMessage(assistantMsg.content);
            return;
          }
        }
        setAssistantMessage(null); // No assistant message found
      } catch (e) {
        setAssistantMessage('Failed to load assistant message.');
        if (e instanceof Error) {
          console.error('Error loading assistant message:', e.message);
        } else {
          console.error('Unknown error loading assistant message:', e);
        }
      } finally {
        setMessageLoading(false);
      }
    };
    fetchAssistantMessage();
  }, [task_id]);

  // Fetch schema when JsonSchema tab is selected
  useEffect(() => {
    if (activeTab === 'JsonSchema' && !toxicitySchema) {
      fetch('/api/schema/toxicity')
        .then(res => res.json())
        .then(setToxicitySchema)
        .catch(() => setToxicitySchema({ error: "Failed to load schema" }));
    }
  }, [activeTab, toxicitySchema]);

  useEffect(() => {
    fetchTask();
  }, [fetchTask]);

  useEffect(() => {
    if (!task_id) return;
    // Fetch files for the task
    fetch(`/api/tasks/${task_id}/files`, { credentials: "include" })
      .then(res => res.json())
      .then(data => setTaskFiles(data.files || []))
      .catch(() => setTaskFiles([]));
  }, [task_id]);

  useEffect(() => {
    if (task && task.environment_id) {
      fetch(`/api/environments/${task.environment_id}/files`, { credentials: "include" })
        .then(res => res.json())
        .then(data => setEnvFiles(data.files || []))
        .catch(() => setEnvFiles([]));
    }
  }, [task]);

  if (loading) return <div className="text-white p-8">Loading...</div>;
  if (!task) return <div className="text-white p-8">Task not found.</div>;

  const duration = getDuration(task.created_at, task.finished_at);

  return (
    <div className="max-w-7xl mx-auto p-8 bg-gray-900 rounded-lg shadow text-white mt-12 flex flex-col min-h-[70vh]">
      <HomeButton className="absolute top-8 left-8" />
      <div className="flex flex-row w-full min-h-[60vh]">
        {/* Sidebar: File List */}
        <div className="w-64 pr-6 border-r border-gray-800 flex-shrink-0 overflow-y-auto" style={{ minWidth: 220 }}>
          <div className="mb-6">
            <div className="mb-2 font-semibold text-green-300">Output files for this Task:</div>
            {taskFiles.length === 0 ? (
              <div className="text-gray-400 mb-2">No files for this task.</div>
            ) : (
              <ul className="mb-2">
                {(() => {
                  const maxFilenameLength = 30;
                  const truncate = (name: string) => {
                    if (name.length <= maxFilenameLength) return name;
                    const startLen = maxFilenameLength - 10;
                    const extMatch = name.match(/(\.[^./\\]+)$/);
                    const ext = extMatch ? extMatch[1] : '';
                    return `${name.slice(0, startLen)}...${ext}`;
                  };
                  return taskFiles.map(f => (
                    <li key={f.file_id} className="flex items-center space-x-2">
                      <button
                        className="text-blue-400 underline hover:text-blue-300 cursor-pointer bg-transparent border-none p-0 m-0"
                        title={f.filename}
                        style={{
                          background: 'none',
                          fontSize: '12px',
                          lineHeight: '1',
                          padding: '0',
                          border: 'none',
                          height: '20px',
                          minHeight: 'unset',
                        }}
                        onClick={() => {
                          setPreviewFileId(f.file_id);
                          setPreviewOpen(true);
                        }}
                      >
                        {truncate(f.filename)}
                      </button>
                      <a
                        href={`/api/tasks/${f.task_id}/files/${f.file_id}/download`}
                        className="ml-2 px-1 py-1 text-gray-300 bg-transparent hover:bg-green-400/10 hover:border-green-300 hover:text-green-200 rounded-full text-xs font-semibold flex items-center gap-1 transition focus:outline-none focus:ring-2 focus:ring-green-400/50"
                        title={`Download ${f.filename}`}
                        download
                      >
                        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 20" fill="currentColor" className="w-5 h-5 mr-1"><path d="M9 3.75a1 1 0 0 1 2 0v6.19l1.72-1.72a1 1 0 1 1 1.42 1.42l-3.43 3.43a1 1 0 0 1-1.42 0l-3.43-3.43a1 1 0 1 1 1.42-1.42L9 9.94V3.75ZM4.25 15a.75.75 0 0 1 .75-.75h10a.75.75 0 0 1 0 1.5H5a.75.75 0 0 1-.75-.75Z" /></svg>
                      </a>
                    </li>
                  ));
                })()}
              </ul>
            )}
            <div className="mb-2 font-semibold text-blue-300">Available files for this Environment:</div>
            {envFiles.length === 0 ? (
              <div className="text-gray-400">No files for this environment.</div>
            ) : (
              <ul>
                {(() => {
                  const maxFilenameLength = 30;
                  const truncate = (name: string) => {
                    if (name.length <= maxFilenameLength) return name;
                    const startLen = maxFilenameLength - 10;
                    const extMatch = name.match(/(\.[^./\\]+)$/);
                    const ext = extMatch ? extMatch[1] : '';
                    return `${name.slice(0, startLen)}...${ext}`;
                  };
                  return envFiles.map(f => (
                    <li key={f.file_id} className="flex items-center space-x-2">
                      <button
                        className="text-blue-400 underline hover:text-blue-300 cursor-pointer bg-transparent border-none p-0 m-0"
                        title={f.filename}
                        style={{
                          background: 'none',
                          fontSize: '12px',
                          lineHeight: '1',
                          padding: '0',
                          border: 'none',
                          height: '20px',
                          minHeight: 'unset',
                        }}
                        onClick={() => {
                          setPreviewFileId(f.file_id);
                          setPreviewOpen(true);
                        }}
                      >
                        {truncate(f.filename)}
                      </button>
                      <a
                        href={`/api/environments/${f.environment_id}/files/${f.file_id}/download`}
                        className="ml-2 px-0 py-1 text-gray-300 bg-transparent hover:bg-green-400/10 hover:border-green-300 hover:text-green-200 rounded-full text-xs font-semibold flex items-center gap-1 transition focus:outline-none focus:ring-2 focus:ring-green-400/50"
                        title={`Download ${f.filename}`}
                        download
                      >
                        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 20" fill="currentColor" className="w-5 h-5 mr-1"><path d="M9 3.75a1 1 0 0 1 2 0v6.19l1.72-1.72a1 1 0 1 1 1.42 1.42l-3.43 3.43a1 1 0 0 1-1.42 0l-3.43-3.43a1 1 0 1 1 1.42-1.42L9 9.94V3.75ZM4.25 15a.75.75 0 0 1 .75-.75h10a.75.75 0 0 1 0 1.5H5a.75.75 0 0 1-.75-.75Z" /></svg>
                      </a>
                    </li>
                  ));
                })()}
              </ul>
            )}
          </div>
        </div>
        {/* Main Content */}
        <div className="flex-1 flex flex-col items-center justify-center w-full pl-8">
          {/* Tab content */}
          <div className="flex-1 flex flex-col items-center justify-center w-full">
            {activeTab === 'MDrender' && (
              <div className="w-full max-w-full bg-gray-800 rounded-lg p-6 shadow-inner">
                {messageLoading ? (
                  <div className="text-gray-400">Loading result...</div>
                ) : assistantMessage ? (
                  <div
                    className="prose prose-invert w-full max-w-full"
                    style={{ lineHeight: 2, maxHeight: 650, overflowY: 'auto' }}
                  >
                    <ReactMarkdown remarkPlugins={[remarkGfm]}>{assistantMessage}</ReactMarkdown>
                  </div>
                ) : (
                  <div className="text-gray-400">No assistant message found for this task.</div>
                )}
              </div>
            )}
            {/* {activeTab === 'MDraw' && (
              <div className="w-full max-w-full bg-gray-800 rounded-lg p-6 shadow-inner">
                <div className="text-lg font-semibold mb-2">Raw Markdown</div>
                <pre
                  style={{
                    whiteSpace: 'pre-wrap',
                    wordBreak: 'break-all',
                    textAlign: 'left',
                    color: '#a3e635',
                    background: 'transparent',
                    fontSize: '1.1em',
                    maxHeight: 650,
                    overflowY: 'auto'
                  }}
                >
                  {assistantMessage}
                </pre>
              </div>
            )}
            {activeTab === 'JsonSchema' && (
              <div
                className="w-full max-w-full bg-gray-800 rounded-lg p-6 shadow-inner"
              >
                <div className="text-lg font-semibold mb-2">Toxicity Schema</div>
                <pre style={{
                  textAlign: 'left',
                  color: '#a3e635',
                  background: 'transparent',
                  fontSize: '1.1em',
                  maxHeight: 650,
                  overflowY: 'auto'
                }}>
                  {toxicitySchema ? JSON.stringify(toxicitySchema, null, 2) : "Loading..."}
                </pre>
              </div>
            )} */}
          </div>
        </div>
      </div>
      <footer className="mt-6 pt-1 pb-1 pl-20 pr-10 border-t border-gray-700 bg-gray-950 rounded-b-lg shadow-inner">
        <div className="w-full max-w-10xl mx-auto">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-x-8 gap-y-0 text-base overflow-x-auto min-w-0 w-full">
            <div className="flex items-center space-x-2">
              <span className="font-semibold text-gray-300">Tool:</span>
              <span>{task.workflow_id === 1 ? "ToxIndex RAP" : task.workflow_id === 2 ? "ToxIndex Vanilla" : task.workflow_id ?? <span className="text-gray-400">Unknown</span>}</span>
            </div>
            <div className="flex items-center space-x-2 min-w-0 w-full">
              <span className="font-semibold text-gray-300 whitespace-nowrap">Executed by user:</span>
              <span className="flex-1 min-w-0 w-full">
                {task.user_id ? (
                  <button
                    className="text-purple-400 underline hover:text-purple-300 focus:outline-none bg-transparent p-0 shadow-none border-none break-all w-full text-left"
                    style={{
                      background: 'none',
                      boxShadow: 'none',
                      border: 'none',
                      wordBreak: 'break-all',
                      minHeight: '1rem',
                      lineHeight: '1rem',
                      paddingTop: '1px',
                      paddingBottom: '1px',
                      display: 'inline-block',
                    }}
                    onClick={() => navigate(`/user/${task.user_id}`)}
                    title="Go to user profile"
                    aria-label="Go to user profile"
                  >
                    {task.user_id}
                  </button>
                ) : (
                  <span className="text-gray-400">Unknown</span>
                )}
              </span>
            </div>
            <div className="flex items-center space-x-2">
              <span className="font-semibold text-gray-300">Created:</span>
              <span>{task.created_at ? new Date(task.created_at).toLocaleString() : <span className="text-gray-400">Unknown</span>}</span>
            </div>
            <div className="flex items-center space-x-2">
              <span className="font-semibold text-gray-300">Finished:</span>
              <span>{task.finished_at ? new Date(task.finished_at).toLocaleString() : <span className="text-gray-400">-</span>}</span>
            </div>
            <div className="flex items-center space-x-2">
              <span className="font-semibold text-gray-300">Duration:</span>
              <span>{duration ?? <span className="text-gray-400">-</span>}</span>
            </div>
            <div className="flex items-center space-x-2 min-w-0 w-full">
              <span className="font-semibold text-gray-300 whitespace-nowrap">Environment:</span>
              <span className="flex-1 min-w-0 w-full">
                {task.environment_id ? (
                  <button
                    className="text-blue-400 underline hover:text-blue-300 focus:outline-none bg-transparent p-0 shadow-none border-none break-all w-full text-left"
                    style={{
                      background: 'none',
                      boxShadow: 'none',
                      border: 'none',
                      wordBreak: 'break-all',
                      minHeight: '1rem',
                      lineHeight: '1rem',
                      paddingTop: '1px',
                      paddingBottom: '1px',
                      display: 'inline-block',
                    }}
                    onClick={() => navigate(`/environment/${task.environment_id}`)}
                    title="Go to environment"
                    aria-label="Go to environment"
                  >
                    {task.environment_id}
                  </button>
                ) : (
                  <span className="text-gray-400">None</span>
                )}
              </span>
            </div>
            <div className="flex items-center space-x-2 min-w-0 w-full">
              <span className="font-semibold text-gray-300 whitespace-nowrap">Chat Session:</span>
              <span className="flex-1 min-w-0 w-full">
                {task.session_id ? (
                  <button
                    className="text-green-400 underline hover:text-green-300 focus:outline-none bg-transparent p-0 shadow-none border-none break-all w-full text-left"
                    style={{
                      background: 'none',
                      boxShadow: 'none',
                      border: 'none',
                      wordBreak: 'break-all',
                      minHeight: '1rem',
                      lineHeight: '1rem',
                      paddingTop: '1px',
                      paddingBottom: '1px',
                      display: 'inline-block',
                    }}
                    onClick={() => navigate(`/chat/session/${task.session_id}`)}
                    title="Go to chat session"
                    aria-label="Go to chat session"
                  >
                    {task.session_id}
                  </button>
                ) : (
                  <span className="text-gray-400">None</span>
                )}
              </span>
            </div>
          </div>
        </div>
      </footer>
      <FilePreviewModal
        fileId={previewFileId}
        isOpen={previewOpen}
        onRequestClose={() => setPreviewOpen(false)}
      />
    </div>
  );
};

export default TaskDetail; 