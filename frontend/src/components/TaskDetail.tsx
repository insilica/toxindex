import React, { useEffect, useState, useCallback } from "react";
import { useParams, useNavigate } from "react-router-dom";
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import HomeButton from "./shared/HomeButton";

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
  const [activeTab, setActiveTab] = useState<'MDrender' | 'MDraw' | 'JsonSchema'>('MDrender');
  const [toxicitySchema, setToxicitySchema] = useState<any>(null);

  const [taskFiles, setTaskFiles] = useState<any[]>([]);
  const [envFiles, setEnvFiles] = useState<any[]>([]);

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
      <div className="mb-6 w-full">
        <div className="mb-2 font-semibold text-green-300">Files for this Task:</div>
        {taskFiles.length === 0 ? (
          <div className="text-gray-400 mb-2">No files for this task.</div>
        ) : (
          <ul className="mb-2">
            {taskFiles.map(f => (
              <li key={f.file_id}>
                <a
                  href={`/api/environments/${f.environment_id}/files/${f.file_id}/download`}
                  className="text-blue-400 underline hover:text-blue-300"
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  {f.filename}
                </a>
              </li>
            ))}
          </ul>
        )}
        <div className="mb-2 font-semibold text-blue-300">Files for this Environment:</div>
        {envFiles.length === 0 ? (
          <div className="text-gray-400">No files for this environment.</div>
        ) : (
          <ul>
            {envFiles.map(f => (
              <li key={f.file_id}>
                <a
                  href={`/api/environments/${f.environment_id}/files/${f.file_id}/download`}
                  className="text-blue-400 underline hover:text-blue-300"
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  {f.filename}
                </a>
              </li>
            ))}
          </ul>
        )}
      </div>
      <div className="flex items-center justify-center mb-5" style={{ minHeight: '2.0rem' }}>
        <h2 className="text-3xl font-semibold text-center mr-6" style={{ marginBottom: 0 }}>{task.title}</h2>
        <div className="flex space-x-2 ml-4">
          <button
            className={`px-4 py-2 font-semibold rounded-t-lg focus:outline-none transition border-b-4 ${activeTab === 'MDrender' ? 'border-green-400 text-green-300 bg-gray-800' : 'border-transparent text-gray-400 bg-gray-900 hover:text-white'}`}
            onClick={() => setActiveTab('MDrender')}
          >
            MDrender
          </button>
          <button
            className={`px-4 py-2 font-semibold rounded-t-lg focus:outline-none transition border-b-4 ${activeTab === 'MDraw' ? 'border-green-400 text-green-300 bg-gray-800' : 'border-transparent text-gray-400 bg-gray-900 hover:text-white'}`}
            onClick={() => setActiveTab('MDraw')}
          >
            MDraw
          </button>
          <button
            className={`px-4 py-2 font-semibold rounded-t-lg focus:outline-none transition border-b-4 ${activeTab === 'JsonSchema' ? 'border-green-400 text-green-300 bg-gray-800' : 'border-transparent text-gray-400 bg-gray-900 hover:text-white'}`}
            onClick={() => setActiveTab('JsonSchema')}
          >
            JsonSchema
          </button>
        </div>
      </div>
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
        {activeTab === 'MDraw' && (
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
        )}
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
    </div>
  );
};

export default TaskDetail; 